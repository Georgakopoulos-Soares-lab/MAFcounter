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
#include <cstdio>       // for snprintf

// ------------------- For 128-bit support ---------------------
#include <boost/multiprecision/cpp_int.hpp>
namespace bmp = boost::multiprecision;
using uint128_t = bmp::uint128_t;
// ------------------------------------------------------------

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

    uint128_t mask = (((uint128_t)1 << (2*K)) - 1);
    uint128_t kmerVal = combined & mask;
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
    // ... (unchanged single-file output code) ...
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
            for(uint8_t i=0;i<nIDs;i++){
                uint16_t len=0;
                if (!m_ifs.read((char*)&len,sizeof(len))) { m_eof=true; return; }
                m_ifs.ignore(len);
            }
        }
        {
            char e[4];
            if (!m_ifs.read(e,4)) { m_eof=true; return; }
            if(std::strncmp(e,"ENDG",4)!=0) { m_eof=true; return; }
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
                uint128_t remainder = (((uint128_t)val[1]) << 56) | ((val[0] >> 8) & 0x00FFFFFFFFFFFFFFULL);
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

// Global variables used in multiple-files mode.
static std::vector<std::string> g_genomeNames; 
static int g_k=0;
static int g_numThreads=1;
static uint8_t g_nIDs=0;

// For partitioning work among threads.
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

// ----- New helper functions for multiple-files mmap output -----
//
// For each output line we write: "<kmer> <count>\n"
// So the line length is: kmer.size() + 1 + digitCount(count) + 1.
static inline size_t computeMultipleLineSize(const std::string &kmer, uint32_t count) {
    return kmer.size() + 1 + digitCount(count) + 1;
}

// Build the output line into 'dest' starting at offset 'written' and update 'written'
static void buildMultipleLine(char *dest,
                              const std::string &kmer,
                              uint32_t count,
                              size_t &written)
{
    std::memcpy(dest + written, kmer.data(), kmer.size());
    written += kmer.size();
    dest[written++] = ' ';
    char cbuf[16];
    int len = std::snprintf(cbuf, sizeof(cbuf), "%u", count);
    std::memcpy(dest + written, cbuf, len);
    written += len;
    dest[written++] = '\n';
}

// ----- Phase 1 worker: compute required output bytes per genome -----
//
// Each thread processes its p1 range and, for every aggregated kmer record,
// for each genome with a nonzero count, adds the line length to its local tally.
static void phase1MultipleMmapWorker(int threadID,
                                     uint16_t p1Start,
                                     uint16_t p1End,
                                     const std::string &finalBinPath,
                                     std::vector<size_t> &localGenomeBytes,
                                     int k,
                                     uint8_t nIDs)
{
    FinalBinReader reader(finalBinPath, k, p1Start, p1End);
    if (!reader.good()) return;
    bool haveCurr = false;
    uint128_t currKmer = 0;
    std::vector<uint32_t> counts(nIDs, 0);
    while (true) {
        NextRecord nr = reader.getNext();
        if (!nr.valid) {
            if (haveCurr) {
                KmerResult r = decodeKmer128(currKmer, k);
                for (uint8_t gid = 0; gid < nIDs; gid++) {
                    if (counts[gid] > 0) {
                        size_t lineSize = computeMultipleLineSize(r.kmer, counts[gid]);
                        localGenomeBytes[gid] += lineSize;
                    }
                }
            }
            break;
        }
        if (!haveCurr) {
            haveCurr = true;
            currKmer = nr.kmer;
            std::fill(counts.begin(), counts.end(), 0);
            counts[nr.genomeID] += nr.count;
        } else {
            if (nr.kmer == currKmer) {
                counts[nr.genomeID] += nr.count;
            } else {
                KmerResult r = decodeKmer128(currKmer, k);
                for (uint8_t gid = 0; gid < nIDs; gid++) {
                    if (counts[gid] > 0) {
                        size_t lineSize = computeMultipleLineSize(r.kmer, counts[gid]);
                        localGenomeBytes[gid] += lineSize;
                    }
                }
                currKmer = nr.kmer;
                std::fill(counts.begin(), counts.end(), 0);
                counts[nr.genomeID] += nr.count;
            }
        }
    }
}

// ----- Phase 2 worker: write output lines into mmapped files -----
//
// Each thread reprocesses its p1 range; for every aggregated record it writes,
// for each genome with nonzero count, the output line into the mmap region at its assigned offset.
static void phase2MultipleMmapWorker(int threadID,
                                     uint16_t p1Start,
                                     uint16_t p1End,
                                     const std::string &finalBinPath,
                                     const std::vector<std::vector<size_t>> &threadGenomeOffsets,
                                     const std::vector<char*> &mappedPtrs,
                                     int k,
                                     uint8_t nIDs)
{
    // Make a local copy of the per-genome starting offsets for this thread.
    std::vector<size_t> offsets = threadGenomeOffsets[threadID];
    FinalBinReader reader(finalBinPath, k, p1Start, p1End);
    if (!reader.good()) return;
    bool haveCurr = false;
    uint128_t currKmer = 0;
    std::vector<uint32_t> counts(nIDs, 0);
    while (true) {
        NextRecord nr = reader.getNext();
        if (!nr.valid) {
            if (haveCurr) {
                KmerResult r = decodeKmer128(currKmer, k);
                for (uint8_t gid = 0; gid < nIDs; gid++) {
                    if (counts[gid] > 0) {
                        buildMultipleLine(mappedPtrs[gid], r.kmer, counts[gid], offsets[gid]);
                    }
                }
            }
            break;
        }
        if (!haveCurr) {
            haveCurr = true;
            currKmer = nr.kmer;
            std::fill(counts.begin(), counts.end(), 0);
            counts[nr.genomeID] += nr.count;
        } else {
            if (nr.kmer == currKmer) {
                counts[nr.genomeID] += nr.count;
            } else {
                KmerResult r = decodeKmer128(currKmer, k);
                for (uint8_t gid = 0; gid < nIDs; gid++) {
                    if (counts[gid] > 0) {
                        buildMultipleLine(mappedPtrs[gid], r.kmer, counts[gid], offsets[gid]);
                    }
                }
                currKmer = nr.kmer;
                std::fill(counts.begin(), counts.end(), 0);
                counts[nr.genomeID] += nr.count;
            }
        }
    }
}

// ----- New doMultipleFilesOutputWithMmap -----
//
// This function first launches a set of threads (each processing its own p1 range)
// to compute the total number of bytes that will be needed per genome (first pass).
// Then it sums the per-thread contributions and allocates one mmapped output file per genome.
// In a second pass the threads re-read their assigned records and write their lines
// into the proper file at the precomputed offsets.
int doMultipleFilesOutputWithMmap(const std::string &finalBinPath, int numThreads)
{
    using clock = std::chrono::high_resolution_clock;
    // First, read header information from finalBinPath to initialize g_k, g_nIDs, and g_genomeNames.
    {
        std::ifstream ifs(finalBinPath, std::ios::binary);
        if (!ifs) {
            std::cerr << "Cannot open " << finalBinPath << "\n";
            return 1;
        }
        char mg[4];
        if (!ifs.read(mg, 4)) { return 1; }
        if (std::strncmp(mg, "BINF", 4) != 0) { return 1; }
        if (!ifs.read((char*)&g_k, sizeof(g_k))) { return 1; }
        if (!ifs.read((char*)&g_nIDs, sizeof(g_nIDs))) { return 1; }
        g_genomeNames.resize(g_nIDs);
        for (uint8_t i = 0; i < g_nIDs; i++) {
            uint16_t len;
            if (!ifs.read((char*)&len, sizeof(len))) { return 1; }
            g_genomeNames[i].resize(len);
            if (!ifs.read((char*)g_genomeNames[i].data(), len)) { return 1; }
        }
        char endg[4];
        if (!ifs.read(endg, 4)) { return 1; }
        if (std::strncmp(endg, "ENDG", 4) != 0) { return 1; }
    }
    
    // Partition p1 values among threads.
    int totalP1 = (1 << P1_BITS);
    if (numThreads < 1) numThreads = 1;
    int chunkSize = (totalP1 + numThreads - 1) / numThreads;
    std::vector<ThreadRange> thrRanges(numThreads);
    for (int i = 0; i < numThreads; i++) {
        uint16_t start = (uint16_t)(i * chunkSize);
        uint16_t end = (uint16_t)std::min<int>((i + 1) * chunkSize - 1, totalP1 - 1);
        thrRanges[i] = {start, end};
    }
    
    // Phase 1: each thread computes local byte counts per genome.
    // Create a per-thread vector of size g_nIDs, initialized to 0.
    std::vector<std::vector<size_t>> localGenomeBytes(numThreads, std::vector<size_t>(g_nIDs, 0));
    {
        std::vector<std::thread> threads;
        threads.reserve(numThreads);
        for (int t = 0; t < numThreads; t++) {
            threads.emplace_back(phase1MultipleMmapWorker,
                                 t,
                                 thrRanges[t].p1Start,
                                 thrRanges[t].p1End,
                                 finalBinPath,
                                 std::ref(localGenomeBytes[t]),
                                 g_k,
                                 g_nIDs);
        }
        for (auto &th : threads)
            th.join();
    }
    
    // Sum per-thread counts for each genome and compute per-thread starting offsets.
    std::vector<size_t> genomeTotalBytes(g_nIDs, 0);
    // threadGenomeOffsets[threadID][genomeID]
    std::vector<std::vector<size_t>> threadGenomeOffsets(numThreads, std::vector<size_t>(g_nIDs, 0));
    for (uint8_t gid = 0; gid < g_nIDs; gid++) {
        size_t offset = 0;
        for (int t = 0; t < numThreads; t++) {
            threadGenomeOffsets[t][gid] = offset;
            offset += localGenomeBytes[t][gid];
        }
        genomeTotalBytes[gid] = offset;
    }
    
    // For each genome, create the output file and mmapped region.
    std::vector<char*> mappedPtrs(g_nIDs, nullptr);
    std::vector<int> fds(g_nIDs, -1);
    for (uint8_t gid = 0; gid < g_nIDs; gid++) {
        std::string outName = "genome_" + g_genomeNames[gid] + ".txt";
        int fd = open(outName.c_str(), O_RDWR | O_CREAT, 0666);
        if (fd < 0) {
            std::cerr << "Cannot create " << outName << "\n";
            return 1;
        }
        if (ftruncate(fd, genomeTotalBytes[gid]) != 0) {
            std::cerr << "ftruncate failed for " << outName << ": " << errno << "\n";
            close(fd);
            return 1;
        }
        void* mapPtr = mmap(nullptr, genomeTotalBytes[gid], PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
        if (mapPtr == MAP_FAILED) {
            std::cerr << "mmap failed for " << outName << ": " << errno << "\n";
            close(fd);
            return 1;
        }
        mappedPtrs[gid] = (char*)mapPtr;
        fds[gid] = fd;
    }
    
    // Phase 2: re-read the records and write them into the mmapped files at the proper offsets.
    {
        std::vector<std::thread> threads;
        threads.reserve(numThreads);
        for (int t = 0; t < numThreads; t++) {
            threads.emplace_back(phase2MultipleMmapWorker,
                                 t,
                                 thrRanges[t].p1Start,
                                 thrRanges[t].p1End,
                                 finalBinPath,
                                 std::cref(threadGenomeOffsets),
                                 std::cref(mappedPtrs),
                                 g_k,
                                 g_nIDs);
        }
        for (auto &th : threads)
            th.join();
    }
    
    // Sync and unmap each genome file.
    for (uint8_t gid = 0; gid < g_nIDs; gid++) {
        if (msync(mappedPtrs[gid], genomeTotalBytes[gid], MS_SYNC) != 0) {
            std::cerr << "msync failed for genome " << gid << "\n";
        }
        munmap(mappedPtrs[gid], genomeTotalBytes[gid]);
        close(fds[gid]);
    }
    
    std::cout << "Created per-genome output using mmap. Byte sizes per genome:\n";
    for (uint8_t gid = 0; gid < g_nIDs; gid++) {
        std::cout << "  " << g_genomeNames[gid] << ": " << genomeTotalBytes[gid] << " bytes\n";
    }
    return 0;
}

// ----------------------
// The rest of the code (e.g., doSingleFileOutputWithMmap) remains unchanged.
// ----------------------

void printUsage() {
    std::cout <<
        "Usage:\n"
        "  ./maf_counter_dump [options] <final.bin>\n\n"
        "Options:\n"
        "  --output_mode=<single|multiple>\n"
        "        Output mode: single = single-file output using mmap (default),\n"
        "                     multiple = per-genome output using mmap (updated version).\n"
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
            finalBin = arg;
        }
    }

    if (finalBin.empty()) {
        std::cerr << "Error: Missing required <final.bin> argument.\n";
        printUsage();
        return 1;
    }

    if (outputMode == "multiple") {
        return doMultipleFilesOutputWithMmap(finalBin, numThreads);
    } else if (outputMode == "single") {
        return doSingleFileOutput(finalBin, outputFile);
    } else {
        std::cerr << "Error: Unknown output mode '" << outputMode << "'.\n";
        printUsage();
        return 1;
    }
}