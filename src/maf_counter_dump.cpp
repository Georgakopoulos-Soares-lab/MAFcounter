#include <iostream>    // for std::cin, std::cout, etc.
#include <fstream>     // for std::ifstream, std::ofstream
#include <cstdint>
#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <sstream>
#include <cstring>
#include <array>
#include <filesystem>
#include <thread>
#include <mutex>
#include <atomic>
#include <chrono>
#include <sstream>
#include <sys/mman.h>   // for mmap
#include <fcntl.h>      // for open
#include <unistd.h>     // for close, ftruncate
#include <errno.h>

// ------------------- For 128-bit support ---------------------
#include <boost/multiprecision/cpp_int.hpp>
namespace bmp = boost::multiprecision;
using uint128_t = bmp::uint128_t;
// ------------------------------------------------------------

// We now allow k up to 64
static const int KMER_MAX_LENGTH = 64;
static const int P1_BITS         = 10;  // used in multiple-files approach
static const int P2_BITS         = 10;  // used in multiple-files approach
static const int MAX_GENOME_IDS  = 256; // upper bound

static const char BASE_LUT[4] = {'A','C','G','T'};

static std::array<std::atomic<size_t>, 1024> g_p1SizesAtomic;
static bool g_decodeByteInited = false;
static std::array<std::string, 256> g_decodeByte;

static void initDecodeByte() {
    for (int x = 0; x < 256; x++) {
        int c0 = (x >> 6) & 3;
        int c1 = (x >> 4) & 3;
        int c2 = (x >> 2) & 3;
        int c3 = x & 3;
        g_decodeByte[x] = {
            BASE_LUT[c0],
            BASE_LUT[c1],
            BASE_LUT[c2],
            BASE_LUT[c3]
        };
    }
}

struct KmerResult {
    std::string kmer;
    uint32_t    count;
};

static KmerResult decodeKmer128(uint128_t combined, int K)
{
    if (!g_decodeByteInited) {
        initDecodeByte();
        g_decodeByteInited = true;
    }

    // The lower 2*K bits are the kmer. The rest (above 2*K bits) could be count or other info
    // but in decoding just the kmer string, we focus on the lower 2*K bits as the kmer.
    uint128_t mask = (((uint128_t)1 << (2*K)) - 1);
    uint128_t kmerVal = combined & mask;
    // For printing or final usage, we might interpret the rest as count in some contexts.
    // But decodeKmer128 is used primarily to decode the base sequence from the lower 2*K bits.
    uint32_t cnt = (uint32_t)(combined >> (2*K)).convert_to<uint64_t>();

    int totalBits = 2 * K;
    int leftoverBits = totalBits % 8;
    int fullBytes = totalBits / 8;

    std::string result;
    result.reserve(K);

    if (leftoverBits > 0) {
        int shift = totalBits - leftoverBits;
        uint64_t partialVal = ((kmerVal >> shift) & (((uint128_t)1 << leftoverBits) - 1)).convert_to<uint64_t>();
        int leftoverBases = leftoverBits / 2;
        for (int i = 0; i < leftoverBases; i++) {
            int baseShift = 2 * (leftoverBases - 1 - i);
            int baseCode = (partialVal >> baseShift) & 3;
            result.push_back(BASE_LUT[baseCode]);
        }
    }
    for (int i = fullBytes - 1; i >= 0; i--) {
        uint64_t shift = i * 8;
        uint64_t byteVal = ((kmerVal >> shift) & 0xFF).convert_to<uint64_t>();
        result += g_decodeByte[(uint8_t)byteVal];
    }
    if ((int)result.size() > K) {
        result.resize(K);
    }
    return { result, cnt };
}

struct FileHeader {
    char     magic[4];
    uint32_t k;
    uint16_t p1_bits;
    uint16_t p2_bits;
    uint8_t  genome_id_bits;
    uint8_t  suffix_mode;  // 0 -> 32-bit suffix, 1 -> 64-bit suffix, 3 -> 128-bit suffix
    uint8_t  reserved[2];
};

int doSingleFileOutput(const std::string &finalBinPath, const std::string &optionalOutName)
{
    std::ifstream ifs(finalBinPath, std::ios::binary);
    if (!ifs) {
        std::cerr << "Cannot open " << finalBinPath << "\n";
        return 1;
    }
    char magic[5];
    magic[4] = 0;
    if (!ifs.read(magic, 4)) {
        std::cerr << "Error reading magic\n";
        return 1;
    }
    if (std::strncmp(magic, "BINF", 4) != 0) {
        std::cerr << "Not a valid binary format\n";
        return 1;
    }
    int k = 0;
    if (!ifs.read((char*)&k, sizeof(k))) {
        std::cerr << "Error reading k\n";
        return 1;
    }
    if (k > KMER_MAX_LENGTH) {
        std::cerr << "Error: k in file is too large\n";
        return 1;
    }
    uint8_t nIDs = 0;
    if (!ifs.read((char*)&nIDs, sizeof(nIDs))) {
        std::cerr << "Error reading nIDs\n";
        return 1;
    }
    std::vector<std::string> genomeNames(nIDs);
    for (uint8_t i = 0; i < nIDs; i++) {
        uint16_t len;
        if (!ifs.read((char*)&len, sizeof(len))) {
            std::cerr << "Error reading genome name length\n";
            return 1;
        }
        genomeNames[i].resize(len);
        if (!ifs.read((char*)genomeNames[i].data(), len)) {
            std::cerr << "Error reading genome name\n";
            return 1;
        }
    }
    {
        char sentinelG[5];
        sentinelG[4] = 0;
        if (!ifs.read(sentinelG, 4)) {
            std::cerr << "Error reading sentinel\n";
            return 1;
        }
        if (std::strncmp(sentinelG, "ENDG", 4) != 0) {
            std::cerr << "Missing ENDG\n";
            return 1;
        }
    }
    std::string outName = optionalOutName;
    if (outName.empty()) {
        outName = "final_sorted_" + std::to_string(k) + "_dump.txt";
    }
    std::ofstream ofs(outName);
    if (!ofs) {
        std::cerr << "Cannot create " << outName << "\n";
        return 1;
    }

    if (k <= 10) {
        // small K mode
        uint64_t mapSize = 0;
        if(!ifs.read((char*)&mapSize, sizeof(mapSize))) {
            std::cerr << "Error reading map size\n";
            return 1;
        }
        uint64_t totalRead = 0;
        while (totalRead < mapSize) {
            uint32_t key32 = 0;
            if (!ifs.read((char*)&key32, sizeof(key32))) {
                break;
            }
            std::array<uint32_t, 256> arr{};
            if (!ifs.read((char*)arr.data(), 256*sizeof(uint32_t))) {
                break;
            }
            totalRead++;
            uint128_t fullVal128 = key32;
            KmerResult r = decodeKmer128(fullVal128, k);
            ofs << r.kmer << ": ";
            bool firstItem = true;
            for (int g = 0; g < 256; g++) {
                if (arr[g] > 0) {
                    if (!firstItem) ofs << ", ";
                    firstItem = false;
                    if (g < (int)genomeNames.size())
                        ofs << genomeNames[g] << ": " << arr[g];
                    else
                        ofs << "UnknownID#" << g << ": " << arr[g];
                }
            }
            ofs << "\n";
        }
    }
    else {
        // big K mode
        while (true) {
            char binMagic[5];
            binMagic[4] = 0;
            if (!ifs.read(binMagic, 4)) {
                break;
            }
            if (std::strncmp(binMagic, "BINc", 4) != 0) {
                break;
            }
            uint16_t p1 = 0;
            if(!ifs.read((char*)&p1, sizeof(p1))) {
                break;
            }
            char kmhdMagic[4];
            if(!ifs.read(kmhdMagic, 4)) {
                break;
            }
            if(std::strncmp(kmhdMagic, "KMHD", 4) != 0) {
                break;
            }
            FileHeader hdr;
            std::memcpy(hdr.magic, kmhdMagic, 4);
            if(!ifs.read((char*)&hdr.k,           sizeof(hdr.k)))           break;
            if(!ifs.read((char*)&hdr.p1_bits,     sizeof(hdr.p1_bits)))     break;
            if(!ifs.read((char*)&hdr.p2_bits,     sizeof(hdr.p2_bits)))     break;
            if(!ifs.read((char*)&hdr.genome_id_bits, sizeof(hdr.genome_id_bits))) break;
            if(!ifs.read((char*)&hdr.suffix_mode, sizeof(hdr.suffix_mode))) break;
            if(!ifs.read((char*)&hdr.reserved[0], 2)) break;

            char sentinel1[4];
            if (!ifs.read(sentinel1, 4)) {
                break;
            }
            uint32_t groupCount = 0;
            if(!ifs.read((char*)&groupCount, sizeof(groupCount))) {
                break;
            }
            struct GMeta {
                uint16_t p2_val;
                uint32_t sz;
            };
            std::vector<GMeta> metas(groupCount);
            for(uint32_t i = 0; i < groupCount; i++){
                if(!ifs.read((char*)&metas[i].p2_val, sizeof(metas[i].p2_val))) break;
                if(!ifs.read((char*)&metas[i].sz,     sizeof(metas[i].sz)))     break;
            }
            char sentinel2[4];
            if(!ifs.read(sentinel2,4)) {
                break;
            }

            int remainder_bits = 2*hdr.k - hdr.p1_bits - hdr.p2_bits;

            auto flushKmer = [&](uint128_t fullVal,
                                 const std::unordered_map<uint8_t,uint32_t>& gcounts)
            {
                KmerResult r = decodeKmer128(fullVal, k);
                ofs << r.kmer << ": ";
                bool first = true;
                for (auto &gc : gcounts){
                    if(!first) ofs << ", ";
                    first = false;
                    if(gc.first < genomeNames.size())
                        ofs << genomeNames[gc.first] << ": " << gc.second;
                    else
                        ofs << "UnknownID#" << (int)gc.first << ": " << gc.second;
                }
                ofs << "\n";
            };

            for(auto &gm : metas){
                if(gm.sz == 0) continue;

                if(hdr.suffix_mode == 0) {
                    // suffix in 32 bits
                    std::vector<uint32_t> arr(gm.sz);
                    if(!ifs.read((char*)arr.data(), gm.sz*sizeof(uint32_t))) {
                        break;
                    }
                    bool first = true;
                    uint128_t currVal = 0;
                    std::unordered_map<uint8_t,uint32_t> gcounts;
                    for(auto v : arr){
                        uint8_t gID   = (uint8_t)(v & 0xFF);
                        uint64_t rem  = (v >> 8);
                        uint128_t fullVal = p1;
                        fullVal <<= (hdr.p2_bits + remainder_bits);
                        fullVal |= ((uint128_t)gm.p2_val << remainder_bits);
                        fullVal |= rem;
                        if(first){
                            first = false;
                            currVal = fullVal;
                            gcounts.clear();
                            gcounts[gID] = 1;
                        } else {
                            if(fullVal == currVal){
                                gcounts[gID]++;
                            } else {
                                flushKmer(currVal, gcounts);
                                currVal = fullVal;
                                gcounts.clear();
                                gcounts[gID] = 1;
                            }
                        }
                    }
                    if(!first){
                        flushKmer(currVal, gcounts);
                    }
                }
                else if(hdr.suffix_mode == 1) {
                    // suffix in 64 bits
                    std::vector<uint64_t> arr(gm.sz);
                    if(!ifs.read((char*)arr.data(), gm.sz*sizeof(uint64_t))) {
                        break;
                    }
                    bool first = true;
                    uint128_t currVal = 0;
                    std::unordered_map<uint8_t,uint32_t> gcounts;
                    for(auto v : arr){
                        uint8_t gID = (uint8_t)(v & 0xFF);
                        uint64_t rem = (v >> 8);
                        uint128_t fullVal = p1;
                        fullVal <<= (hdr.p2_bits + remainder_bits);
                        fullVal |= ((uint128_t)gm.p2_val << remainder_bits);
                        fullVal |= rem;
                        if(first){
                            first = false;
                            currVal = fullVal;
                            gcounts.clear();
                            gcounts[gID] = 1;
                        } else {
                            if(fullVal == currVal){
                                gcounts[gID]++;
                            } else {
                                flushKmer(currVal, gcounts);
                                currVal = fullVal;
                                gcounts.clear();
                                gcounts[gID] = 1;
                            }
                        }
                    }
                    if(!first){
                        flushKmer(currVal, gcounts);
                    }
                }
                else if(hdr.suffix_mode == 3) {
                    // suffix in 128 bits
                    // We'll read each record as two uint64_t = total 128 bits.
                    std::vector<std::array<uint64_t,2>> arr(gm.sz);
                    if(!ifs.read((char*)arr.data(), gm.sz*2*sizeof(uint64_t))) {
                        break;
                    }
                    bool first = true;
                    uint128_t currVal = 0;
                    std::unordered_map<uint8_t,uint32_t> gcounts;
                    for(auto &v : arr){
                        // lower 8 bits of v[0] = genome ID
                        uint8_t gID = (uint8_t)(v[0] & 0xFF);
                        // next 120 bits are remainder
                        uint128_t rem = (( (uint128_t)v[1]) << 56) | ((v[0] >> 8) & 0x00FFFFFFFFFFFFFFULL);

                        uint128_t fullVal = p1;
                        fullVal <<= (hdr.p2_bits + remainder_bits);
                        fullVal |= ((uint128_t)gm.p2_val << remainder_bits);
                        fullVal |= rem;

                        if(first) {
                            first = false;
                            currVal = fullVal;
                            gcounts.clear();
                            gcounts[gID] = 1;
                        } else {
                            if(fullVal == currVal) {
                                gcounts[gID]++;
                            } else {
                                flushKmer(currVal, gcounts);
                                currVal = fullVal;
                                gcounts.clear();
                                gcounts[gID] = 1;
                            }
                        }
                    }
                    if(!first){
                        flushKmer(currVal, gcounts);
                    }
                }
                else {
                    // unknown suffix_mode
                    std::cerr << "Unknown suffix_mode encountered.\n";
                    return 1;
                }
            }
            char bend[5];
            bend[4] = 0;
            if(!ifs.read(bend,4)){
                break;
            }
            if(std::strncmp(bend,"BEND",4)!=0){
                break;
            }
        }
    }
    ofs.close();
    return 0;
}

struct NextRecord {
    bool      valid;
    uint128_t kmer;
    uint8_t   genomeID;
    uint32_t  count;
};

static inline uint128_t reassembleKmer(uint16_t p1, uint16_t p2, uint128_t remainder, int k)
{
    int remainder_bits = 2*k - P1_BITS - P2_BITS;
    uint128_t fullVal = p1;
    fullVal <<= (P2_BITS + remainder_bits);
    fullVal |= ((uint128_t)p2 << remainder_bits);
    fullVal |= remainder;
    return fullVal;
}

static inline uint16_t getP1Range(uint128_t fullKmer, int k)
{
    int totalBits = 2*k;
    int shift = totalBits - P1_BITS;
    if (shift < 0) return 0;
    uint128_t mask = (((uint128_t)1 << P1_BITS) - 1);
    uint128_t p1val = (fullKmer >> shift) & mask;
    return static_cast<uint16_t>(p1val.convert_to<uint64_t>());
}

class FinalBinReader {
public:
    FinalBinReader(const std::string &fname, int k, uint16_t p1Start, uint16_t p1End)
        : m_k(k), m_p1Start(p1Start), m_p1End(p1End)
    {
        m_ifs.open(fname, std::ios::binary);
        if (!m_ifs) {
            m_eof = true;
            return;
        }
        {
            char mg[4];
            if (!m_ifs.read(mg,4)) { m_eof = true; return; }
            if (std::strncmp(mg,"BINF",4)!=0) {
                m_eof=true; return;
            }
        }
        {
            int fileK=0;
            if (!m_ifs.read((char*)&fileK,sizeof(fileK))) { m_eof=true; return; }
        }
        {
            uint8_t nIDs=0;
            if (!m_ifs.read((char*)&nIDs,sizeof(nIDs))) { m_eof=true; return; }
            for(uint8_t i=0; i<nIDs; i++){
                uint16_t len=0;
                if (!m_ifs.read((char*)&len,sizeof(len))) { m_eof=true; return; }
                m_ifs.ignore(len);
            }
        }
        {
            char e[4];
            if (!m_ifs.read(e,4)) { m_eof=true; return; }
            if(std::strncmp(e,"ENDG",4)!=0) {
                m_eof=true; return;
            }
        }

        m_smallKmode = (k <= 10);
        if (m_smallKmode) {
            if (!m_ifs.read((char*)&m_mapSize,sizeof(m_mapSize))) {
                m_eof=true; 
                return;
            }
            m_mapReadCount = 0;
            m_idx=256;
        } else {
            // We assume there's a "final.metadata" for big K
            std::ifstream meta("final.metadata");
            if(!meta) {
                m_eof=true;
                return;
            }
            std::vector<std::pair<uint16_t,uint64_t>> metaEntries;
            std::string line;
            while(std::getline(meta,line)){
                if(line.empty()) continue;
                std::istringstream iss(line);
                uint16_t mp1; 
                uint64_t off;
                if(!(iss>>mp1>>off)) continue;
                metaEntries.push_back({mp1,off});
            }
            meta.close();
            std::sort(metaEntries.begin(), metaEntries.end(), 
                      [](auto&a,auto&b){return a.first<b.first;});
            auto it = std::lower_bound(metaEntries.begin(), metaEntries.end(), m_p1Start,
                [](const std::pair<uint16_t,uint64_t>& e, uint16_t val){
                    return e.first<val;
                }
            );
            if(it==metaEntries.end()|| it->first>m_p1End){
                m_eof=true;
                return;
            }
            uint64_t seekOffset = it->second;
            m_ifs.seekg(seekOffset);
            if(!m_ifs){
                m_eof=true;
                return;
            }
            m_inBinBlock = false;
        }
    }

    bool good() const { return !m_eof; }

    NextRecord getNext()
    {
        if(m_eof) return {false,0,0,0};
        if(m_smallKmode){
            return getNextSmall();
        } else {
            return getNextBig();
        }
    }

private:
    NextRecord getNextSmall()
    {
        while(true){
            if(m_eof) return {false,0,0,0};
            if(m_mapReadCount>=m_mapSize){
                char endm[4];
                m_ifs.read(endm,4);
                m_eof=true;
                return {false,0,0,0};
            }
            if(m_idx>=256){
                uint32_t key32=0;
                if(!m_ifs.read((char*)&key32,sizeof(key32))){
                    m_eof=true;
                    return {false,0,0,0};
                }
                if(!m_ifs.read((char*)m_counts.data(),256*sizeof(uint32_t))){
                    m_eof=true;
                    return {false,0,0,0};
                }
                m_mapReadCount++;
                m_idx=0;
                m_currKmer = key32;
            }
            while(m_idx<256 && m_counts[m_idx]==0){
                m_idx++;
            }
            if(m_idx>=256) continue;
            NextRecord nr;
            nr.valid = true;
            nr.kmer = m_currKmer;
            nr.genomeID = (uint8_t)m_idx;
            nr.count = m_counts[m_idx];
            m_idx++;
            uint16_t p1 = getP1Range(nr.kmer,m_k);
            if(p1>=m_p1Start && p1<=m_p1End){
                return nr;
            }
        }
    }

    NextRecord getNextBig()
    {
        while(true){
            if(m_eof)return{false,0,0,0};
            if(!m_inBinBlock){
                char binC[4];
                if(!m_ifs.read(binC,4)){ m_eof=true; return {false,0,0,0};}
                if(std::strncmp(binC,"BINc",4)!=0){m_eof=true;return{false,0,0,0};}
                if(!m_ifs.read((char*)&m_blockP1,sizeof(m_blockP1))){m_eof=true;return{false,0,0,0};}
                if(m_blockP1>m_p1End){m_eof=true;return{false,0,0,0};}
                char kmhdMagic[4];
                if(!m_ifs.read(kmhdMagic,4)){m_eof=true;return{false,0,0,0};}
                if(std::strncmp(kmhdMagic,"KMHD",4)!=0){m_eof=true;return{false,0,0,0};}
                FileHeader hdr;
                std::memcpy(hdr.magic, kmhdMagic,4);
                if(!m_ifs.read((char*)&hdr.k,sizeof(hdr.k))){m_eof=true;return{false,0,0,0};}
                if(!m_ifs.read((char*)&hdr.p1_bits,sizeof(hdr.p1_bits))){m_eof=true;return{false,0,0,0};}
                if(!m_ifs.read((char*)&hdr.p2_bits,sizeof(hdr.p2_bits))){m_eof=true;return{false,0,0,0};}
                if(!m_ifs.read((char*)&hdr.genome_id_bits,sizeof(hdr.genome_id_bits))){m_eof=true;return{false,0,0,0};}
                if(!m_ifs.read((char*)&hdr.suffix_mode,sizeof(hdr.suffix_mode))){m_eof=true;return{false,0,0,0};}
                if(!m_ifs.read((char*)&hdr.reserved[0],2)){m_eof=true;return{false,0,0,0};}
                char s1[4];
                if(!m_ifs.read(s1,4)){m_eof=true;return{false,0,0,0};}
                uint32_t groupCount=0;
                if(!m_ifs.read((char*)&groupCount,sizeof(groupCount))){m_eof=true;return{false,0,0,0};}
                m_groups.clear();
                m_groups.resize(groupCount);
                for(uint32_t i=0;i<groupCount;i++){
                    if(!m_ifs.read((char*)&m_groups[i].p2,sizeof(m_groups[i].p2))){m_eof=true;return{false,0,0,0};}
                    if(!m_ifs.read((char*)&m_groups[i].sz,sizeof(m_groups[i].sz))){m_eof=true;return{false,0,0,0};}
                }
                char s2[4];
                if(!m_ifs.read(s2,4)){m_eof=true;return{false,0,0,0};}
                m_suffixMode = hdr.suffix_mode;
                m_inBinBlock=true;
            }
            while(m_groupIndex<m_groups.size() && m_groups[m_groupIndex].done()){
                m_groupIndex++;
            }
            if(m_groupIndex>=m_groups.size()){
                char bend[4];
                if(!m_ifs.read(bend,4)){m_eof=true;return{false,0,0,0};}
                m_inBinBlock=false;
                m_groupIndex=0;
                continue;
            }
            auto &cg = m_groups[m_groupIndex];
            if(!cg.loaded){
                if(m_suffixMode==0){
                    cg.data32.resize(cg.sz);
                    if(!m_ifs.read((char*)cg.data32.data(),cg.sz*sizeof(uint32_t))) {m_eof=true;return{false,0,0,0};}
                } else if(m_suffixMode==1) {
                    cg.data64.resize(cg.sz);
                    if(!m_ifs.read((char*)cg.data64.data(),cg.sz*sizeof(uint64_t))) {m_eof=true;return{false,0,0,0};}
                } else if(m_suffixMode==3) {
                    cg.data128.resize(cg.sz);
                    if(!m_ifs.read((char*)cg.data128.data(), cg.sz*2*sizeof(uint64_t))) {m_eof=true;return{false,0,0,0};}
                } else {
                    m_eof=true;
                    return {false,0,0,0};
                }
                cg.loaded=true;
                cg.idx=0;
            }
            if(m_suffixMode==0){
                if(cg.idx>=cg.data32.size()){cg.idx=cg.sz;continue;}
                uint32_t val = cg.data32[cg.idx++];
                uint8_t gID=(uint8_t)(val & 0xFF);
                uint64_t remainder=(val>>8);
                uint128_t fullKmer = reassembleKmer(m_blockP1,cg.p2,remainder,m_k);
                uint16_t p1 = getP1Range(fullKmer,m_k);
                if(p1>=m_p1Start && p1<=m_p1End){
                    return {true,fullKmer,gID,1};
                }
            } else if(m_suffixMode==1){
                if(cg.idx>=cg.data64.size()){cg.idx=cg.sz;continue;}
                uint64_t v=cg.data64[cg.idx++];
                uint8_t gID=(uint8_t)(v&0xFF);
                uint64_t remainder=(v>>8);
                uint128_t fullKmer = reassembleKmer(m_blockP1,cg.p2,remainder,m_k);
                uint16_t p1=getP1Range(fullKmer,m_k);
                if(p1>=m_p1Start && p1<=m_p1End){
                    return {true,fullKmer,gID,1};
                }
            } else if(m_suffixMode==3){
                if(cg.idx>=cg.data128.size()){cg.idx=cg.sz;continue;}
                auto &val = cg.data128[cg.idx++];
                uint8_t gID = (uint8_t)(val[0] & 0xFF);
                uint128_t remainder = (( (uint128_t)val[1]) << 56) | ((val[0] >> 8) & 0x00FFFFFFFFFFFFFFULL);
                uint128_t fullKmer = reassembleKmer(m_blockP1, cg.p2, remainder, m_k);
                uint16_t p1=getP1Range(fullKmer,m_k);
                if(p1>=m_p1Start && p1<=m_p1End){
                    return {true,fullKmer,gID,1};
                }
            }
        }
    }

    std::ifstream m_ifs;
    bool m_eof=false;
    bool m_smallKmode=false;

    uint64_t m_mapSize=0;
    uint64_t m_mapReadCount=0;
    std::array<uint32_t,256> m_counts{};
    uint128_t m_currKmer=0;
    int m_idx=256;
    uint8_t m_suffixMode=0;

    bool m_inBinBlock=false;
    uint16_t m_blockP1=0;
    struct GroupData {
        uint16_t p2=0;
        uint32_t sz=0;
        bool loaded=false;
        std::vector<uint64_t> data64;
        std::vector<uint32_t> data32;
        std::vector<std::array<uint64_t,2>> data128;
        size_t idx=0;
        bool done() const {return idx>=sz;}
    };
    std::vector<GroupData> m_groups;
    size_t m_groupIndex=0;

    int m_k;
    uint16_t m_p1Start;
    uint16_t m_p1End;
};

static std::vector<std::string> g_genomeNames; 
static int g_k=0;
static int g_numThreads=1;
static uint8_t g_nIDs=0;

struct ThreadRange {
    uint16_t p1Start;
    uint16_t p1End;
};

static inline int digitCount(uint32_t num) {
    if (num < 10) return 1;
    if (num < 100) return 2;
    if (num < 1000) return 3;
    if (num < 10000) return 4;
    if (num < 100000) return 5;
    if (num < 1000000) return 6;
    if (num < 10000000) return 7;
    if (num < 100000000) return 8;
    if (num < 1000000000) return 9;
    return 10;  
}

static size_t computeLineSize(const std::string &kmerStr,
                              const std::unordered_map<uint8_t, uint32_t> &gcounts,
                              const std::vector<std::string> &genomeNames)
{
    size_t lineSize = kmerStr.size() + 2; 
    bool first = true;
    for(auto &kv : gcounts){
        if(kv.second==0) continue;
        if(!first) lineSize+=2;
        first=false;
        uint8_t gID=kv.first;
        if(gID<genomeNames.size()) {
            lineSize+= genomeNames[gID].size();
        } else {
            lineSize+= 9;
            if(gID>=100) lineSize+=3;
            else if(gID>=10) lineSize+=2;
            else lineSize+=1;
        }
        lineSize += 2;
        lineSize += digitCount(kv.second);
    }
    lineSize += 1;
    return lineSize;
}

static void p1SizeComputerWorker(int threadID,
                                 uint16_t p1Start,
                                 uint16_t p1End,
                                 const std::string &finalBinPath,
                                 const std::vector<std::string> &genomeNames,
                                 int k,
                                 uint8_t nIDs)
{
    FinalBinReader reader(finalBinPath, k, p1Start, p1End);
    if(!reader.good()) return;
    bool haveCurr=false;
    uint128_t currKmer=0;
    std::vector<uint32_t> counts(nIDs,0);
    auto flushKmer = [&](uint128_t km){
        std::unordered_map<uint8_t,uint32_t> gmap;
        for(uint8_t gid=0;gid<nIDs;gid++){
            if(counts[gid]>0){
                gmap[gid]=counts[gid];
            }
        }
        if(gmap.empty()) return;
        KmerResult r=decodeKmer128(km,k);
        size_t lineBytes=computeLineSize(r.kmer,gmap,genomeNames);
        uint16_t p1Val=getP1Range(km,k);
        g_p1SizesAtomic[p1Val].fetch_add(lineBytes,std::memory_order_relaxed);
    };
    while(true){
        NextRecord nr=reader.getNext();
        if(!nr.valid){
            if(haveCurr) flushKmer(currKmer);
            break;
        }
        if(!haveCurr){
            haveCurr=true;
            currKmer=nr.kmer;
            std::fill(counts.begin(),counts.end(),0);
            counts[nr.genomeID]+=nr.count;
        } else {
            if(nr.kmer==currKmer){
                counts[nr.genomeID]+=nr.count;
            } else {
                flushKmer(currKmer);
                currKmer=nr.kmer;
                std::fill(counts.begin(),counts.end(),0);
                counts[nr.genomeID]+=nr.count;
            }
        }
    }
}

int computeLineSizesSingleFileKmer(const std::string &finalBinPath, int numThreads)
{
    std::ifstream ifs(finalBinPath, std::ios::binary);
    if(!ifs){
        std::cerr<<"Cannot open "<<finalBinPath<<"\n";
        return 1;
    }
    {
        char mg[4];
        if(!ifs.read(mg,4)){return 1;}
        if(std::strncmp(mg,"BINF",4)!=0){
            return 1;
        }
    }
    int k=0;
    if(!ifs.read((char*)&k,sizeof(k))){
        return 1;
    }
    if(k<=0 || k>KMER_MAX_LENGTH){
        return 1;
    }
    uint8_t nIDs=0;
    if(!ifs.read((char*)&nIDs,sizeof(nIDs))){
        return 1;
    }
    std::vector<std::string> genomeNames(nIDs);
    for(uint8_t i=0;i<nIDs;i++){
        uint16_t len=0;
        if(!ifs.read((char*)&len,sizeof(len))){return 1;}
        genomeNames[i].resize(len);
        if(!ifs.read((char*)genomeNames[i].data(),len)){return 1;}
    }
    {
        char endg[4];
        if(!ifs.read(endg,4)){return 1;}
        if(std::strncmp(endg,"ENDG",4)!=0){
            return 1;
        }
    }
    ifs.close();
    for(auto &x:g_p1SizesAtomic){
        x.store(0,std::memory_order_relaxed);
    }
    int totalP1=(1<<P1_BITS);
    if(numThreads<1) numThreads=1;
    int chunkSize=(totalP1+numThreads-1)/numThreads;
    std::vector<ThreadRange> thrRanges(numThreads);
    for(int i=0;i<numThreads;i++){
        uint16_t start=(uint16_t)(i*chunkSize);
        uint16_t end=(uint16_t)std::min<int>((i+1)*chunkSize-1,totalP1-1);
        thrRanges[i]={start,end};
    }
    std::vector<std::thread> threads;
    threads.reserve(numThreads);
    for(int t=0;t<numThreads;t++){
        threads.emplace_back(p1SizeComputerWorker,t,
                             thrRanges[t].p1Start,thrRanges[t].p1End,
                             finalBinPath,std::cref(genomeNames),k,nIDs);
    }
    for(auto &th:threads){
        th.join();
    }
    size_t totalBytes=0;
    for(int i=0;i<1024;i++){
        totalBytes+= g_p1SizesAtomic[i].load(std::memory_order_relaxed);
    }
    std::cout<<"[computeLineSizesSingleFileKmer] Estimated total bytes = "<<totalBytes<<"\n";
    return 0;
}

// ----------------------------------------------------------------------------------
// doSingleFileOutputWithMmap (uses line-size computation + parallel write)
// ----------------------------------------------------------------------------------

static std::vector<size_t> g_p1Offsets(1025);

static void buildLine(char *dest,
                      const std::string &kmerStr,
                      const std::unordered_map<uint8_t,uint32_t> &gcounts,
                      const std::vector<std::string> &genomeNames,
                      size_t &written)
{
    // We'll directly build: "kmer: name: count, name2: count2, ...\n"
    // We'll track 'written' as we go.

    // 1) copy kmerStr
    std::memcpy(dest + written, kmerStr.data(), kmerStr.size());
    written += kmerStr.size();
    // 2) ": "
    dest[written++] = ':';
    dest[written++] = ' ';
    // 3) each item
    bool first=true;
    for(auto &kv:gcounts){
        if(kv.second==0) continue;
        if(!first){
            dest[written++] = ',';
            dest[written++] = ' ';
        }
        first=false;
        uint8_t gID=kv.first;
        if(gID<genomeNames.size()){
            // copy genome name
            const std::string &nm= genomeNames[gID];
            std::memcpy(dest+written,nm.data(),nm.size());
            written+= nm.size();
        } else {
            // "UnknownID#xxx"
            static const char unknownPrefix[]="UnknownID#";
            std::memcpy(dest+written, unknownPrefix, 9);
            written+=9;
            char buf[4];
            int len=std::snprintf(buf,4,"%u",(unsigned)gID);
            std::memcpy(dest+written,buf,len);
            written+=len;
        }
        // ": "
        dest[written++]=':';
        dest[written++]=' ';
        // count
        char cbuf[16];
        int len=std::snprintf(cbuf,16,"%u",kv.second);
        std::memcpy(dest+written,cbuf,len);
        written+=len;
    }
    // newline
    dest[written++]='\n';
}

static void p1WriteWorker(int threadID,
                          uint16_t p1Start,
                          uint16_t p1End,
                          const std::string &finalBinPath,
                          const std::vector<std::string> &genomeNames,
                          int k,
                          uint8_t nIDs,
                          char *mappedPtr)
{
    FinalBinReader reader(finalBinPath, k, p1Start, p1End);
    if(!reader.good()) return;

    bool haveCurr=false;
    uint128_t currKmer=0;
    std::vector<uint32_t> counts(nIDs,0);

    auto flushKmer=[&](uint128_t km){
        std::unordered_map<uint8_t,uint32_t> gmap;
        for(uint8_t gid=0;gid<nIDs;gid++){
            if(counts[gid]>0) gmap[gid]=counts[gid];
        }
        if(gmap.empty()) return;
        KmerResult r=decodeKmer128(km,k);
        size_t lineLen=computeLineSize(r.kmer,gmap,genomeNames);
        uint16_t p1Val=getP1Range(km,k);
        size_t offset=g_p1Offsets[p1Val];
        buildLine(mappedPtr, r.kmer, gmap, genomeNames, offset);
        g_p1Offsets[p1Val]=offset;
    };

    while(true){
        NextRecord nr=reader.getNext();
        if(!nr.valid){
            if(haveCurr) flushKmer(currKmer);
            break;
        }
        if(!haveCurr){
            haveCurr=true;
            currKmer=nr.kmer;
            std::fill(counts.begin(),counts.end(),0);
            counts[nr.genomeID]+=nr.count;
        } else {
            if(nr.kmer==currKmer){
                counts[nr.genomeID]+=nr.count;
            } else {
                flushKmer(currKmer);
                currKmer=nr.kmer;
                std::fill(counts.begin(),counts.end(),0);
                counts[nr.genomeID]+=nr.count;
            }
        }
    }
}

int doSingleFileOutputWithMmap(const std::string &finalBinPath,
                               const std::string &optionalOutName,
                               int numThreads)
{
    std::ifstream ifs(finalBinPath, std::ios::binary);
    if(!ifs){
        std::cerr<<"Cannot open "<<finalBinPath<<"\n";
        return 1;
    }
    {
        char mg[4];
        if(!ifs.read(mg,4)) {
            std::cerr<<"Error reading magic\n";
            return 1;
        }
        if(std::strncmp(mg,"BINF",4)!=0){
            std::cerr<<"Not a BINF\n";
            return 1;
        }
    }
    int k=0;
    if(!ifs.read((char*)&k,sizeof(k))){
        std::cerr<<"Error reading k\n";
        return 1;
    }
    if(k> KMER_MAX_LENGTH){
        std::cerr<<"k in file is too large\n";
        return 1;
    }
    uint8_t nIDs=0;
    if(!ifs.read((char*)&nIDs,sizeof(nIDs))){
        std::cerr<<"Error reading nIDs\n";
        return 1;
    }
    std::vector<std::string> genomeNames(nIDs);
    for(uint8_t i=0;i<nIDs;i++){
        uint16_t len=0;
        if(!ifs.read((char*)&len,sizeof(len))){
            std::cerr<<"Error reading genome name length\n";
            return 1;
        }
        genomeNames[i].resize(len);
        if(!ifs.read((char*)genomeNames[i].data(),len)){
            std::cerr<<"Error reading genome name\n";
            return 1;
        }
    }
    {
        char endg[4];
        if(!ifs.read(endg,4)){
            std::cerr<<"Error reading ENDG\n";
            return 1;
        }
        if(std::strncmp(endg,"ENDG",4)!=0){
            std::cerr<<"Missing ENDG\n";
            return 1;
        }
    }
    ifs.close();

    for(auto &x:g_p1SizesAtomic){
        x.store(0,std::memory_order_relaxed);
    }
    int totalP1=(1<<P1_BITS);
    if(numThreads<1) numThreads=1;
    int chunkSize=(totalP1+numThreads-1)/numThreads;
    std::vector<ThreadRange> thrRanges(numThreads);
    for(int i=0;i<numThreads;i++){
        uint16_t start=(uint16_t)(i*chunkSize);
        uint16_t end= (uint16_t)std::min<int>((i+1)*chunkSize-1,totalP1-1);
        thrRanges[i]={start,end};
    }

    // First pass: compute line sizes
    {
        std::vector<std::thread> threads;
        threads.reserve(numThreads);
        for(int t=0;t<numThreads;t++){
            threads.emplace_back(p1SizeComputerWorker,t,
                                 thrRanges[t].p1Start,thrRanges[t].p1End,
                                 finalBinPath,std::cref(genomeNames),
                                 k,nIDs);
        }
        for(auto &th:threads){
            th.join();
        }
    }
    // Compute prefix sums
    size_t totalBytes=0;
    for(int i=0;i<1024;i++){
        size_t sz=g_p1SizesAtomic[i].load(std::memory_order_relaxed);
        g_p1Offsets[i]=totalBytes;
        totalBytes+=sz;
    }
    g_p1Offsets[1024]=totalBytes;

    std::string outName= optionalOutName;
    if(outName.empty()){
        outName="final_sorted_"+std::to_string(k)+"_dump.txt";
    }

    // Create the file
    {
        int fd=open(outName.c_str(), O_RDWR|O_CREAT, 0666);
        if(fd<0){
            std::cerr<<"Cannot create "<<outName<<"\n";
            return 1;
        }
        if(ftruncate(fd, totalBytes)!=0){
            std::cerr<<"ftruncate failed: "<<errno<<"\n";
            close(fd);
            return 1;
        }
        void* mapPtr=mmap(nullptr, totalBytes, PROT_READ|PROT_WRITE,MAP_SHARED,fd,0);
        if(mapPtr==MAP_FAILED){
            std::cerr<<"mmap failed: "<<errno<<"\n";
            close(fd);
            return 1;
        }
        char *mappedPtr=(char*)mapPtr;

        // Second pass: write lines
        {
            std::vector<std::thread> threads;
            threads.reserve(numThreads);
            for(int t=0;t<numThreads;t++){
                threads.emplace_back(p1WriteWorker,t,
                                     thrRanges[t].p1Start,thrRanges[t].p1End,
                                     finalBinPath,std::cref(genomeNames),
                                     k,nIDs,
                                     mappedPtr);
            }
            for(auto &th:threads){
                th.join();
            }
        }

        msync(mappedPtr,totalBytes,MS_SYNC);
        munmap(mappedPtr,totalBytes);
        close(fd);
    }
    std::cout<<"Created single-file text output of size "<<totalBytes<<" bytes.\n";
    return 0;
}

// ------------------------------------------------------------
// doMultipleFilesOutput
// ------------------------------------------------------------

// We remove the old "writeKmerCount" approach that merged (kmer+count) into 128 bits.
// Instead, we store them separately as a struct with 128-bit kmer + 32-bit count.

#pragma pack(push, 1)
struct KmerCountRecord {
    uint64_t kmerLo;
    uint64_t kmerHi;
    uint32_t count;
};
#pragma pack(pop)

static std::vector<std::vector<std::ofstream>> g_outBin;

static void phase1Worker(int threadID, ThreadRange range, const std::string &finalBin)
{
    FinalBinReader reader(finalBin,g_k,range.p1Start, range.p1End);
    if(!reader.good()) return;

    bool haveCurr=false;
    uint128_t currKmer=0;
    std::vector<uint32_t> counts(g_nIDs,0);

    auto flushCurr=[&](uint128_t km){
        for(uint8_t gid=0; gid<g_nIDs; gid++){
            uint32_t c=counts[gid];
            if(c>0){
                KmerCountRecord rec;
                rec.kmerLo = (uint64_t)(km & 0xFFFFFFFFFFFFFFFFULL);
                rec.kmerHi = (uint64_t)(km >> 64);
                rec.count  = c;
                g_outBin[threadID][gid].write((const char*)&rec, sizeof(rec));
            }
        }
    };

    while(true){
        NextRecord nr=reader.getNext();
        if(!nr.valid){
            if(haveCurr) flushCurr(currKmer);
            break;
        }
        if(!haveCurr){
            haveCurr=true;
            currKmer=nr.kmer;
            std::fill(counts.begin(),counts.end(),0);
            counts[nr.genomeID]+=nr.count;
        } else {
            if(nr.kmer==currKmer){
                counts[nr.genomeID]+=nr.count;
            } else {
                flushCurr(currKmer);
                currKmer=nr.kmer;
                std::fill(counts.begin(),counts.end(),0);
                counts[nr.genomeID]+=nr.count;
            }
        }
    }
}

static void phase2Worker(int threadID, int numThreads)
{
    for(int gid=threadID; gid<g_nIDs; gid+=numThreads){
        std::string txtName="genome_"+g_genomeNames[gid]+".txt";
        std::ofstream ofs(txtName);
        if(!ofs) {
            continue;
        }
        for(int t=0;t<g_numThreads;t++){
            std::string binName="genome_"+std::to_string(gid)
                               +"_worker_"+std::to_string(t)+".bin";
            if(!std::filesystem::exists(binName)) continue;
            std::ifstream ifb(binName,std::ios::binary);
            if(!ifb) continue;

            while(true){
                KmerCountRecord rec;
                if(!ifb.read((char*)&rec,sizeof(rec))) break;
                uint128_t kmer = ((uint128_t)rec.kmerHi << 64) | rec.kmerLo;
                KmerResult r=decodeKmer128(kmer,g_k);
                ofs<<r.kmer<<" "<<rec.count<<"\n";
            }
            std::error_code ec;
            std::filesystem::remove(binName,ec);
        }
    }
}

int doMultipleFilesOutput(const std::string &finalBinPath, int numThreads)
{
    using clock=std::chrono::high_resolution_clock;
    auto t0=clock::now();
    g_numThreads=numThreads;
    if(g_numThreads<1) g_numThreads=1;

    {
        std::ifstream ifs(finalBinPath, std::ios::binary);
        if(!ifs){
            std::cerr<<"Cannot open "<<finalBinPath<<"\n";
            return 1;
        }
        {
            char mg[4];
            if(!ifs.read(mg,4)){
                return 1;
            }
            if(std::strncmp(mg,"BINF",4)!=0){
                return 1;
            }
        }
        if(!ifs.read((char*)&g_k,sizeof(g_k))){
            return 1;
        }
        if(!ifs.read((char*)&g_nIDs,sizeof(g_nIDs))){
            return 1;
        }
        g_genomeNames.resize(g_nIDs);
        for(uint8_t i=0;i<g_nIDs;i++){
            uint16_t len;
            if(!ifs.read((char*)&len,sizeof(len))){
                return 1;
            }
            g_genomeNames[i].resize(len);
            if(!ifs.read((char*)g_genomeNames[i].data(), len)){
                return 1;
            }
        }
        {
            char endg[4];
            if(!ifs.read(endg,4)){
                return 1;
            }
            if(std::strncmp(endg,"ENDG",4)!=0){
                return 1;
            }
        }
    }
    auto t00=clock::now();

    g_outBin.resize(g_numThreads);
    for(int t=0;t<g_numThreads;t++){
        g_outBin[t].resize(g_nIDs);
        for(uint8_t gid=0;gid<g_nIDs;gid++){
            std::string fname="genome_"+std::to_string(gid)
                             +"_worker_"+std::to_string(t)+".bin";
            g_outBin[t][gid].open(fname,std::ios::binary);
            if(!g_outBin[t][gid]){
                std::cerr<<"Error creating "<<fname<<"\n";
                return 1;
            }
        }
    }
    auto t000=clock::now();

    int totalP1=(1<<P1_BITS);
    int chunkSize=(totalP1+g_numThreads-1)/g_numThreads;
    std::vector<ThreadRange> thrRanges(g_numThreads);
    for(int i=0;i<g_numThreads;i++){
        uint16_t start=(uint16_t)(i*chunkSize);
        uint16_t end=(uint16_t)std::min<int>((i+1)*chunkSize-1,totalP1-1);
        thrRanges[i]={start,end};
    }

    auto tPhase1Start=clock::now();
    {
        std::vector<std::thread> threads;
        threads.reserve(g_numThreads);
        for(int t=0;t<g_numThreads;t++){
            threads.emplace_back(phase1Worker,t,thrRanges[t],finalBinPath);
        }
        for(auto &th:threads){
            th.join();
        }
    }
    auto tPhase1End=clock::now();
    {
        for(int t=0;t<g_numThreads;t++){
            for(uint8_t gid=0;gid<g_nIDs;gid++){
                g_outBin[t][gid].flush();
            }
        }
    }
    auto tx=clock::now();

    auto tPhase2Start=clock::now();
    {
        std::vector<std::thread> threads;
        threads.reserve(g_numThreads);
        for(int t=0;t<g_numThreads;t++){
            threads.emplace_back(phase2Worker,t,g_numThreads);
        }
        for(auto &th:threads){
            th.join();
        }
    }
    auto tPhase2End=clock::now();

    auto tCleanupStart=clock::now();
    // (any additional cleanup if needed)
    auto tCleanupEnd=clock::now();

    auto phase1sec=std::chrono::duration<double>(tPhase1End-tPhase1Start).count();
    auto phase2sec=std::chrono::duration<double>(tPhase2End-tPhase2Start).count();
    auto cleansec=std::chrono::duration<double>(tCleanupEnd-tCleanupStart).count();
    auto totalSec=std::chrono::duration<double>(tCleanupEnd-t0).count();
    auto debug2=std::chrono::duration<double>(t000-t00).count();
    auto debug=std::chrono::duration<double>(tx-tPhase1End).count();

    std::cout<<"Phase 1 took: "<<phase1sec<<" seconds\n";
    std::cout<<"Phase 2 took: "<<phase2sec<<" seconds\n";
    std::cout<<"Cleanup took: "<<cleansec<<" seconds\n";
    std::cout<<"Total took: "<<totalSec<<" seconds\n";
    std::cout<<"Debug: "<<debug<<" seconds\n";
    std::cout<<"Debug2: "<<debug2<<" seconds\n";
    std::cout<<"Done. Created per-genome .txt outputs.\n";
    return 0;
}

// Print updated usage
void printUsage() {
    std::cout <<
        "Usage:\n"
        "  ./maf_counter_dump [options] <final.bin>\n\n"
        "Options:\n"
        "  --output_mode=<single|multiple>\n"
        "        Output mode: single = single-file output using mmap (default),\n"
        "                     multiple = per-genome output.\n"
        "  --threads <VAL>           Number of threads to use (default: 1).\n"
        "  --output_file <filename>  Output file name for single-file mode (default: final_sorted_<k>_dump.txt).\n"
        "  --output_directory <DIR>  Directory for per-genome output files (default: current directory).\n"
        "  --temp_files_dir <DIR>    Directory for intermediate files (default: current directory).\n"
        "\n"
        "  <final.bin>              Input binary file generated by maf_counter_count.\n";
}

int main(int argc, char* argv[]) {
    // Default parameter values
    std::string outputMode = "single";  // Default mode is single (mmap-based)
    int numThreads = 1;
    std::string outputFile;  // For single-file mode; if empty, will be set later based on k.
    std::string outputDirectory = ".";
    std::string tempFilesDir = ".";
    std::string finalBin;

    // Simple command-line parsing
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg.find("--output_mode=") == 0) {
            outputMode = arg.substr(14);
        } else if (arg == "--threads") {
            if (++i < argc) {
                numThreads = std::stoi(argv[i]);
            }
        } else if (arg == "--output_file") {
            if (++i < argc) {
                outputFile = argv[i];
            }
        } else if (arg == "--output_directory") {
            if (++i < argc) {
                outputDirectory = argv[i];
            }
        } else if (arg == "--temp_files_dir") {
            if (++i < argc) {
                tempFilesDir = argv[i];
            }
        } else if (arg == "--help" || arg == "-h") {
            printUsage();
            return 0;
        } else if (arg[0] == '-') {
            std::cerr << "Unknown option: " << arg << "\n";
            printUsage();
            return 1;
        } else {
            // Assume this is the final.bin file
            finalBin = arg;
        }
    }

    if (finalBin.empty()) {
        std::cerr << "Error: Missing required <final.bin> argument.\n";
        printUsage();
        return 1;
    }

    // Call the appropriate routine based on the output mode
    if (outputMode == "multiple") {
        return doMultipleFilesOutput(finalBin, numThreads);
    } else if (outputMode == "single") {
        return doSingleFileOutputWithMmap(finalBin, outputFile, numThreads);
    } else {
        std::cerr << "Error: Unknown output mode '" << outputMode << "'.\n";
        printUsage();
        return 1;
    }
}