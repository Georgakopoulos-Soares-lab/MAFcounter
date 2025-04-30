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
#include <fcntl.h>      // for open
#include <unistd.h>     // for close, ftruncate
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>

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
// ------------------------------------------------------------------------
// Single-threaded, buffered dump for k <= 10 final.bin format
// ------------------------------------------------------------------------
static const size_t OUT_BUF_SIZE = 4 << 20; // 4 MiB buffer



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

            // ----------------------------------------------------------
            // Have we finished the current 256-counter slab?
            // If yes, either load the next key or declare true EOF.
            // ----------------------------------------------------------
            if (m_idx >= 256) {
                    /* Reached the end of the current counter array.
                     * If we have already read every key, this *really* is EOF;
                     * otherwise read the next <key + 256 counters> block.    */
    
                    if (m_mapReadCount == m_mapSize) {
                        char endm[4];
                        m_ifs.read(endm, 4);          // "ENDM"
                        m_eof = true;
                        return {false,0,0,0};
                    }
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

// ------------------------------------------------------------
// Globals for write-pass coordination
// ------------------------------------------------------------
static std::mutex    g_writeMux;
static size_t        g_writeCursor = 0;

// Helper: build one output line into a std::string
static std::string buildLineString(const std::string &kmer,
                                   const std::unordered_map<uint8_t,uint32_t> &gcounts)
{
    std::string line;
    // reserve roughly: k + 2 + (#entries)*(name+4+digits+2) + 1
    line.reserve(kmer.size() + 2 + gcounts.size()*16 + 1);
    line += kmer;
    line += ": ";
    bool first = true;
    for (auto &kv : gcounts) {
        if (kv.second == 0) continue;
        if (!first) line += ", ";
        first = false;
        uint8_t gid = kv.first;
        if (gid < g_genomeNames.size()) {
            line += g_genomeNames[gid];
        } else {
            line += "UnknownID#";
            line += std::to_string(gid);
        }
        line += ": ";
        line += std::to_string(kv.second);
    }
    line += '\n';
    return line;
}


static std::vector<std::pair<uint16_t,uint64_t>> readP1Metadata(const std::string& metadataFile) {
    std::vector<std::pair<uint16_t,uint64_t>> meta;
    std::ifstream metaIn(metadataFile);
    if (!metaIn) {
        std::cerr << "Error opening metadata file: " << metadataFile << "\n";
        return meta;
    }
    uint16_t p1;
    uint64_t off;
    while (metaIn >> p1 >> off) {
        meta.emplace_back(p1, off);
    }
    std::sort(meta.begin(), meta.end(),
              [](auto &a, auto &b){ return a.first < b.first; });
    return meta;
}

// ------------------------------------------------------------------------
// Worker #1 (prediction pass): for each assigned p1, seek into final.bin
// and compute the byte‐size of every line *without* calling decodeKmer128.
// ------------------------------------------------------------------------
static void p1SizePredWorker(int threadID,
                             const std::vector<uint16_t>& p1List,
                             const std::string& finalBinPath,
                             int k, uint8_t nIDs)
{
    for (uint16_t p1 : p1List) {
        FinalBinReader reader(finalBinPath, k, p1, p1);
        if (!reader.good()) continue;

        bool haveCurr = false;
        uint128_t currKmer = 0;
        std::vector<uint32_t> counts(nIDs, 0);

        // same logic as before, but compute length by arithmetic only
        auto flush = [&](uint128_t) {
            size_t lineBytes = k + 2; // kmer + ": "
            bool first = true;
            for (uint8_t gid = 0; gid < nIDs; ++gid) {
                uint32_t c = counts[gid];
                if (c == 0) continue;
                if (!first) lineBytes += 2; // ", "
                first = false;
                const std::string &nm = (gid < g_genomeNames.size()
                                         ? g_genomeNames[gid]
                                         : "UnknownID#" + std::to_string(gid));
                lineBytes += nm.size() + 2;  // name + ": "
                lineBytes += digitCount(c);
            }
            lineBytes += 1; // newline
            g_p1SizesAtomic[p1].fetch_add(lineBytes,
                                          std::memory_order_relaxed);
        };

        while (true) {
            NextRecord nr = reader.getNext();
            if (!nr.valid) {
                if (haveCurr) flush(currKmer);
                break;
            }
            if (!haveCurr) {
                haveCurr = true;
                currKmer = nr.kmer;
                std::fill(counts.begin(), counts.end(), 0);
                counts[nr.genomeID] += nr.count;
            } else if (nr.kmer == currKmer) {
                counts[nr.genomeID] += nr.count;
            } else {
                flush(currKmer);
                currKmer = nr.kmer;
                std::fill(counts.begin(), counts.end(), 0);
                counts[nr.genomeID] += nr.count;
            }
        }
    }
}

// ------------------------------------------------------------------------
// Revised computeLineSizesSingleFileKmer: 
//   • reads the header (k, nIDs, names… same as before)
//   • reads & sorts final.metadata
//   • round-robins p1 IDs across threads
//   • runs prediction workers
//   • computes g_p1Offsets[] and totalBytes
// ------------------------------------------------------------------------
int computeLineSizesSingleFileKmer(const std::string &finalBinPath,
                                   int numThreads)
{
    // --- parse header (exactly as your old code did) ---
    std::ifstream ifs(finalBinPath, std::ios::binary);
    if (!ifs) {
        std::cerr << "Cannot open " << finalBinPath << "\n";
        return 1;
    }
    // --- parse header so we get k, nIDs and genome names ---
    char magic[4];
    if (!ifs.read(magic,4) || std::strncmp(magic,"BINF",4)!=0) {
        std::cerr<<"Bad magic in "<<finalBinPath<<"\n";
        return 1;
    }
    int k = 0;
    if (!ifs.read((char*)&k,sizeof(k)) || k <= 0 || k > KMER_MAX_LENGTH) {
        std::cerr<<"Bad k in header\n";
        return 1;
    }
    uint8_t nIDs = 0;
    if (!ifs.read((char*)&nIDs,sizeof(nIDs))) {
        std::cerr<<"Bad nIDs in header\n";
        return 1;
    }
    g_genomeNames.resize(nIDs);
    for (uint8_t i = 0; i < nIDs; ++i) {
        uint16_t len = 0;
        if (!ifs.read((char*)&len,sizeof(len))) {
            std::cerr<<"Error reading genome name length\n";
            return 1;
        }
        g_genomeNames[i].resize(len);
        if (!ifs.read(g_genomeNames[i].data(),len)) {
            std::cerr<<"Error reading genome name\n";
            return 1;
        }
    }
    char endg[4];
    if (!ifs.read(endg,4) || std::strncmp(endg,"ENDG",4)!=0) {
        std::cerr<<"Missing ENDG sentinel\n";
        return 1;
    }


    // --- read & sort metadata ---
    auto meta = readP1Metadata("final.metadata");
    if (meta.empty()) {
        std::cerr << "No metadata entries found\n";
        return 1;
    }

    // --- distribute p1 values round-robin to each thread ---
    std::vector<std::vector<uint16_t>> threadP1(numThreads);
    for (size_t i = 0; i < meta.size(); ++i) {
        threadP1[i % numThreads].push_back(meta[i].first);
    }

    // reset atomics
    for (auto &x : g_p1SizesAtomic)
        x.store(0, std::memory_order_relaxed);

    // launch prediction threads
    std::vector<std::thread> threads;
    threads.reserve(numThreads);
    for (int t = 0; t < numThreads; ++t) {
        threads.emplace_back(p1SizePredWorker, t,
                             std::cref(threadP1[t]),
                             std::cref(finalBinPath),
                             k, nIDs);
    }
    for (auto &th : threads) th.join();

    // build offsets & compute total
    size_t totalBytes = 0;
    for (auto &e : meta) {
        uint16_t p1 = e.first;
        g_p1Offsets[p1] = totalBytes;
        totalBytes += g_p1SizesAtomic[p1]
                          .load(std::memory_order_relaxed);
    }
    // sentinel at end:
    g_p1Offsets[meta.size()] = totalBytes;

    std::cout << "[computeLineSizesSingleFileKmer] Estimated total bytes = "
              << totalBytes << "\n";
    return 0;
}

// ------------------------------------------------------------------------
// Worker #2: batch each P1 into one big buffer, then write it in one shot.
// ------------------------------------------------------------------------
static void p1WriteWorker(int threadID,
    const std::vector<uint16_t> &p1List,
    const std::string &finalBinPath,
    int k,
    uint8_t nIDs,
    const std::string &outName)
{
    std::ofstream ofs(outName, std::ios::in|std::ios::out|std::ios::binary);
    if (!ofs) {
        std::cerr<<"Cannot open "<<outName<<" for writing\n";
        return;
    }

    for (uint16_t p1 : p1List) {
        // 1) grab this region’s starting offset
        size_t partSize = g_p1SizesAtomic[p1].load(std::memory_order_relaxed);

        // --- Allocate and reserve exactly the needed bytes for this P1 block ---
        std::string p1Buffer;
        p1Buffer.reserve(partSize);

        size_t writeOffset;
        {
            std::lock_guard<std::mutex> lk(g_writeMux);
            writeOffset   = g_writeCursor;
            g_writeCursor += partSize;
        }

        // 2) stream through only that p1
        FinalBinReader reader(finalBinPath, k, p1, p1);
        if (!reader.good()) continue;

        bool haveCurr = false;
        uint128_t currKmer = 0;
        std::vector<uint32_t> counts(nIDs,0);

        // buildLineString() returns the full line including newline
        auto flush = [&](uint128_t km) {
            // collect non-zero counts
            std::unordered_map<uint8_t,uint32_t> gmap;
            for (uint8_t g = 0; g < nIDs; ++g) {
                if (counts[g] > 0) gmap[g] = counts[g];
            }
            if (gmap.empty()) return;

            // build the single line
            KmerResult r = decodeKmer128(km, k);
            std::string line = buildLineString(r.kmer, gmap);

            // append into our big buffer
            p1Buffer += line;
        };

        // read & flush groups
        while (true) {
            NextRecord nr = reader.getNext();
            if (!nr.valid) {
                if (haveCurr) flush(currKmer);
                break;
            }
            if (!haveCurr) {
                haveCurr = true;
                currKmer = nr.kmer;
                std::fill(counts.begin(), counts.end(), 0);
                counts[nr.genomeID] += nr.count;
            }
            else if (nr.kmer == currKmer) {
                counts[nr.genomeID] += nr.count;
            }
            else {
                flush(currKmer);
                currKmer = nr.kmer;
                std::fill(counts.begin(), counts.end(), 0);
                counts[nr.genomeID] += nr.count;
            }
        }

        // 3) write the entire P1 block in one syscall
        ofs.seekp(writeOffset, std::ios::beg);
        ofs.write(p1Buffer.data(), p1Buffer.size());
    }
}
int doSingleFileSmallKmer(const std::string &finalBinPath,
    const std::string &outName)
{
// 1) Open and parse header
std::ifstream ifs(finalBinPath, std::ios::binary);
if (!ifs) {
std::cerr << "Cannot open " << finalBinPath << "\n";
return 1;
}

// magic
char magic[4];
ifs.read(magic, 4);
if (std::strncmp(magic, "BINF", 4) != 0) {
std::cerr << "Bad magic in " << finalBinPath << "\n";
return 1;
}

// k
int k = 0;
ifs.read(reinterpret_cast<char*>(&k), sizeof(k));

// nIDs
uint8_t nIDs = 0;
ifs.read(reinterpret_cast<char*>(&nIDs), sizeof(nIDs));

// --- READ genome names into g_genomeNames! ---
g_genomeNames.clear();
g_genomeNames.resize(nIDs);
for (uint8_t i = 0; i < nIDs; ++i) {
uint16_t len = 0;
ifs.read(reinterpret_cast<char*>(&len), sizeof(len));
std::string nm(len, '\0');
ifs.read(&nm[0], len);
g_genomeNames[i] = std::move(nm);
}

// skip ENDG sentinel
ifs.ignore(4);

// skip the mapSize (uint64_t)
uint64_t mapSize = 0;
ifs.read(reinterpret_cast<char*>(&mapSize), sizeof(mapSize));

// 2) Prepare reader and output
FinalBinReader reader(finalBinPath, k, 0, (1u << P1_BITS) - 1);
if (!reader.good()) {
std::cerr << "Error initializing reader\n";
return 1;
}

std::ofstream ofs(outName);
if (!ofs) {
std::cerr << "Cannot create " << outName << "\n";
return 1;
}

std::string buffer;
buffer.reserve(OUT_BUF_SIZE);

// 3) Iterate NextRecord, grouping by kmer
uint128_t currKmer = 0;
bool haveCurr = false;
std::vector<uint32_t> counts(nIDs, 0);

auto flush_line = [&]() {
KmerResult kr = decodeKmer128(currKmer, k);
buffer += kr.kmer;
buffer += ": ";
bool first = true;
for (uint8_t gid = 0; gid < nIDs; ++gid) {
if (counts[gid] == 0){
     continue;}
if (!first) buffer += ", ";
first = false;
buffer += g_genomeNames[gid];
buffer += ": ";
buffer += std::to_string(counts[gid]);
}
buffer += '\n';
if (buffer.size() >= OUT_BUF_SIZE) {
ofs << buffer;
buffer.clear();
}
};

while (true) {
NextRecord nr = reader.getNext();
if (!nr.valid) {
if (haveCurr) flush_line();
break;
}
if (!haveCurr || nr.kmer != currKmer) {
if (haveCurr) flush_line();
currKmer = nr.kmer;
std::fill(counts.begin(), counts.end(), 0);
haveCurr = true;
}
counts[nr.genomeID] += nr.count;
}

if (!buffer.empty()) ofs << buffer;
std::cout << "Created single-file dump: " << outName << "\n";
return 0;
}



// ------------------------------------------------------------------------
// Updated doSingleFileOutputWithMmap → single-file write using threads + ofs
// ------------------------------------------------------------------------
int doSingleFileOutputWithMmap(const std::string &finalBinPath,
    const std::string &optionalOutName,
    int numThreads)
{
// --- 1) parse header exactly as before to get k, nIDs, g_genomeNames[] ---
std::ifstream ifs(finalBinPath, std::ios::binary);
if (!ifs) {
std::cerr<<"Cannot open "<<finalBinPath<<"\n";
return 1;
}
    // --- parse header so we get k, nIDs and genome names ---
    char magic[4];
    if (!ifs.read(magic,4) || std::strncmp(magic,"BINF",4)!=0) {
        std::cerr<<"Bad magic in "<<finalBinPath<<"\n";
        return 1;
    }
    int k = 0;
    if (!ifs.read((char*)&k,sizeof(k)) || k <= 0 || k > KMER_MAX_LENGTH) {
        std::cerr<<"Bad k in header\n";
        return 1;
    }
    uint8_t nIDs = 0;
    if (!ifs.read((char*)&nIDs,sizeof(nIDs))) {
        std::cerr<<"Bad nIDs in header\n";
        return 1;
    }
    g_genomeNames.resize(nIDs);
    for (uint8_t i = 0; i < nIDs; ++i) {
        uint16_t len = 0;
        if (!ifs.read((char*)&len,sizeof(len))) {
            std::cerr<<"Error reading genome name length\n";
            return 1;
        }
        g_genomeNames[i].resize(len);
        if (!ifs.read(g_genomeNames[i].data(),len)) {
            std::cerr<<"Error reading genome name\n";
            return 1;
        }
    }
    char endg[4];
    if (!ifs.read(endg,4) || std::strncmp(endg,"ENDG",4)!=0) {
        std::cerr<<"Missing ENDG sentinel\n";
        return 1;
    }


// --- 2) prediction pass (same as Part 1) ---
if (computeLineSizesSingleFileKmer(finalBinPath, numThreads) != 0)
return 1;

// --- 3) decide output file name ---
std::string outName = optionalOutName.empty()
 ? ("final_sorted_" + std::to_string(k) + "_dump.txt")
 : optionalOutName;

// --- 4) create & truncate to total size ---
// total size = sum of all g_p1SizesAtomic[p1] = g_writeCursor after computeLineSizes,
// or recompute as:
size_t totalBytes = 0;
auto meta = readP1Metadata("final.metadata");
for (auto &e : meta)
totalBytes += g_p1SizesAtomic[e.first].load(std::memory_order_relaxed);

int fd = open(outName.c_str(), O_RDWR|O_CREAT, 0666);
if (fd < 0) {
std::cerr<<"Cannot create "<<outName<<"\n";
return 1;
}
if (ftruncate(fd, totalBytes) != 0) {
std::cerr<<"ftruncate failed\n";
close(fd);
return 1;
}
close(fd);

// --- 5) reset global cursor & build per-thread p1 lists (round-robin) ---
g_writeCursor = 0;
std::vector<std::vector<uint16_t>> threadP1(numThreads);
for (size_t i=0; i<meta.size(); ++i)
threadP1[i % numThreads].push_back(meta[i].first);

// --- 6) launch write threads ---
std::vector<std::thread> threads;
threads.reserve(numThreads);
for (int t = 0; t < numThreads; ++t) {
threads.emplace_back(p1WriteWorker, t,
  std::cref(threadP1[t]),
  std::cref(finalBinPath),
  k, nIDs,
  std::cref(outName));
}
for (auto &th : threads) th.join();

std::cout<<"Created single-file text output of size "
<<totalBytes<<" bytes.\n";
return 0;
}




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

// ------------------------------------------------------------------------
// New globals for per-genome multi-file output
static std::array<std::mutex, MAX_GENOME_IDS>             g_fileMutexes;
static std::array<std::atomic<size_t>, MAX_GENOME_IDS>    g_fileSizesAtomic;
static std::array<std::atomic<size_t>, MAX_GENOME_IDS>    g_fileOffsets;

static std::array<int, MAX_GENOME_IDS> g_outFds;







// Print updated usage
void printUsage() {
    std::cout <<
        "Usage:\n"
        "  ./maf_counter_dump [options] <final.bin>\n\n"
        "Options:\n"
        "  --output_mode=<single|multiple>\n"
        "        Output mode: single = single-file  (default),\n"
        "                     multiple = per-genome output.\n"
        "  --threads <VAL>           Number of threads to use (default: 1).\n"
        "  --output_file <filename>  Output file name for single-file mode (default: final_sorted_<k>_dump.txt).\n"
        "  --output_directory <DIR>  Directory for per-genome output files (default: current directory).\n"
        "  --temp_files_dir <DIR>    Directory for intermediate files (default: current directory).\n"
        "\n"
        "  <final.bin>              Input binary file generated by maf_counter_count.\n";
}

int doMultipleFilesOutput(const std::string &finalBinPath,
    int numThreads,
    const std::string &outputDirectory)
{
// --- 1) parse header (same as before) ---
std::ifstream ifs(finalBinPath, std::ios::binary);
if (!ifs) {
std::cerr << "Cannot open " << finalBinPath << "\n";
return 1;
}
char mg[4];
if (!ifs.read(mg,4) || std::strncmp(mg,"BINF",4)!=0) {
std::cerr<<"Bad magic in "<<finalBinPath<<"\n";
return 1;
}
if (!ifs.read((char*)&g_k,sizeof(g_k))) return 1;
if (!ifs.read((char*)&g_nIDs,sizeof(g_nIDs))) return 1;
g_genomeNames.resize(g_nIDs);
for(uint8_t i=0;i<g_nIDs;i++){
uint16_t len;
if(!ifs.read((char*)&len,sizeof(len))) return 1;
g_genomeNames[i].resize(len);
if(!ifs.read(g_genomeNames[i].data(), len)) return 1;
}
char endg[4];
if(!ifs.read(endg,4) || std::strncmp(endg,"ENDG",4)!=0){
std::cerr<<"Missing ENDG sentinel\n"; return 1;
}

// --- 2) read & sort metadata & distribute P1s ---
auto meta = readP1Metadata("final.metadata");
if (meta.empty()) {
std::cerr << "No metadata entries found\n";
return 1;
}
// round-robin P1 → threads
std::vector<std::vector<uint16_t>> thrRanges(numThreads);
for (size_t i = 0; i < meta.size(); ++i) {
thrRanges[i % numThreads].push_back(meta[i].first);
}

for (uint8_t gid = 0; gid < g_nIDs; ++gid) {
    g_fileSizesAtomic[gid].store(0, std::memory_order_relaxed);
    g_fileOffsets  [gid].store(0, std::memory_order_relaxed);
}


auto estimateWorker = [&](const std::vector<uint16_t> &p1List){
for (uint16_t p1 : p1List) {
FinalBinReader reader(finalBinPath, g_k, p1, p1);
if (!reader.good()) continue;
bool have = false;
uint128_t curr=0;
std::vector<uint32_t> counts(g_nIDs,0);
auto flushCounts = [&](){
// compute bytes for each genomeID with count>0
size_t base = g_k + 1; // kmer + " "
for (uint8_t gid=0; gid<g_nIDs; ++gid) {
    uint32_t c = counts[gid];
    if (!c) continue;
    // line is: "<kmer> <count>\n"
    size_t sz = base
              + digitCount(c)
              + 1; // '\n'
    g_fileSizesAtomic[gid].fetch_add(sz,
      std::memory_order_relaxed);
}

};
while (true) {
NextRecord nr = reader.getNext();
if (!nr.valid) {
if (have) flushCounts();
break;
}
if (!have) {
have = true;
curr = nr.kmer;
std::fill(counts.begin(), counts.end(), 0);
}
if (nr.kmer == curr) {
counts[nr.genomeID] += nr.count;
} else {
flushCounts();
curr = nr.kmer;
std::fill(counts.begin(), counts.end(), 0);
counts[nr.genomeID] += nr.count;
}
}
}
};

// launch estimation threads
std::vector<std::thread> estThreads;
for (int t = 0; t < numThreads; ++t) {
estThreads.emplace_back(estimateWorker, std::cref(thrRanges[t]));
}
for (auto &th : estThreads) th.join();

// --- 4) TRUNCATE each genome file to its exact size ---

for (uint8_t gid = 0; gid < g_nIDs; ++gid) {
    std::string fn = outputDirectory + "/genome_" + g_genomeNames[gid] + ".txt";
    int fd = open(fn.c_str(), O_RDWR|O_CREAT, 0666);
    if (fd < 0) { /* error */ }
    size_t needed = g_fileSizesAtomic[gid].load();
    if (ftruncate(fd, needed) != 0) { /* error */ }
    g_outFds[gid] = fd;
}


// --- 5) SECOND PASS: buffered write into each genome file ---
const size_t BUF_THRESHOLD = 5 * 1024 * 1024;  // 20 MiB


auto writeWorker = [&](int threadID){
auto &plist = thrRanges[threadID];
// per-genome in-memory buffers
std::vector<std::string> buffers(g_nIDs);
for (uint16_t p1 : plist) {
FinalBinReader reader(finalBinPath, g_k, p1, p1);
if (!reader.good()) continue;
bool have = false;
uint128_t curr=0;
std::vector<uint32_t> counts(g_nIDs,0);
auto flushLine = [&](){
    KmerResult kr = decodeKmer128(curr, g_k);
    for (uint8_t gid=0; gid<g_nIDs; ++gid) {
        uint32_t c = counts[gid];
        if (!c) continue;
        // build "kmer count\n"
        auto &buf = buffers[gid];
        buf += kr.kmer;
        buf += ' ';
        buf += std::to_string(c);
        buf += '\n';
    }
};
;
;
while (true) {
NextRecord nr = reader.getNext();
if (!nr.valid) {
if (have) flushLine();
break;
}
if (!have) {
have = true;
curr = nr.kmer;
std::fill(counts.begin(), counts.end(), 0);
}
if (nr.kmer == curr) {
counts[nr.genomeID] += nr.count;
} else {
flushLine();
curr = nr.kmer;
std::fill(counts.begin(), counts.end(), 0);
counts[nr.genomeID] += nr.count;
}
uint8_t gid = nr.genomeID;
auto &b = buffers[gid];
if (b.size() >= BUF_THRESHOLD) {
  std::lock_guard<std::mutex> lk(g_fileMutexes[gid]);
  size_t off = g_fileOffsets[gid].fetch_add(b.size());
  pwrite(g_outFds[gid], b.data(), b.size(), off);
  b.clear();
}
}
}
// final flush of any leftovers
    for (uint8_t gid=0; gid<g_nIDs; ++gid) {
          auto &b = buffers[gid];
          if (b.empty()) continue;
          std::lock_guard<std::mutex> lk(g_fileMutexes[gid]);
          size_t off = g_fileOffsets[gid].fetch_add(b.size());
          pwrite(g_outFds[gid], b.data(), b.size(), off);
          b.clear();
        }
};

// launch writer threads
std::vector<std::thread> wThreads;
for (int t = 0; t < numThreads; ++t) {
wThreads.emplace_back(writeWorker, t);
}
for (auto &th : wThreads) th.join();


for (uint8_t gid = 0; gid < g_nIDs; ++gid) {
    close(g_outFds[gid]);
}

return 0;
}

int doMultipleFileSmallKmer(const std::string &finalBinPath,
    const std::string &outputDirectory)
{
// 1) parse header
std::ifstream ifs(finalBinPath, std::ios::binary);
if (!ifs) {
std::cerr << "Cannot open " << finalBinPath << "\n";
return 1;
}
char magic[4];
if (!ifs.read(magic,4) || std::strncmp(magic,"BINF",4)!=0) {
std::cerr << "Bad magic in " << finalBinPath << "\n";
return 1;
}
int k = 0;
if (!ifs.read(reinterpret_cast<char*>(&k), sizeof(k))) {
std::cerr << "Error reading k\n";
return 1;
}
uint8_t nIDs = 0;
if (!ifs.read(reinterpret_cast<char*>(&nIDs), sizeof(nIDs))) {
std::cerr << "Error reading nIDs\n";
return 1;
}
g_genomeNames.resize(nIDs);
for (uint8_t i = 0; i < nIDs; ++i) {
uint16_t len = 0;
if (!ifs.read(reinterpret_cast<char*>(&len), sizeof(len))) {
std::cerr << "Error reading genome name length\n";
return 1;
}
g_genomeNames[i].resize(len);
if (!ifs.read(g_genomeNames[i].data(), len)) {
std::cerr << "Error reading genome name\n";
return 1;
}
}
char endg[4];
if (!ifs.read(endg,4) || std::strncmp(endg,"ENDG",4)!=0) {
std::cerr << "Missing ENDG sentinel\n";
return 1;
}
// skip the mapSize for small‐k format
uint64_t mapSize = 0;
if (!ifs.read(reinterpret_cast<char*>(&mapSize), sizeof(mapSize))) {
std::cerr << "Error reading mapSize\n";
return 1;
}

// 2) open one output file per genome and keep them open
std::vector<std::ofstream> ofs(nIDs);
for (uint8_t gid = 0; gid < nIDs; ++gid) {
std::string fn = outputDirectory + "/genome_" + g_genomeNames[gid] + ".txt";
ofs[gid].open(fn, std::ios::binary);
if (!ofs[gid]) {
std::cerr << "Cannot create " << fn << "\n";
return 1;
}
}

// 3) per‐genome buffers of 4 MiB
std::vector<std::string> buffers(nIDs);
for (auto &b : buffers) b.reserve(OUT_BUF_SIZE);

// 4) single‐threaded scan & write
FinalBinReader reader(finalBinPath, k, 0, (1u << P1_BITS) - 1);
if (!reader.good()) {
std::cerr << "Error initializing reader\n";
return 1;
}

bool haveCurr = false;
uint128_t currKmer = 0;
std::vector<uint32_t> counts(nIDs, 0);

auto flushLine = [&]() {
KmerResult kr = decodeKmer128(currKmer, k);
for (uint8_t gid = 0; gid < nIDs; ++gid) {
uint32_t c = counts[gid];
if (c == 0) continue;
auto &buf = buffers[gid];
buf += kr.kmer;
buf += ' ';
buf += std::to_string(c);
buf += '\n';
if (buf.size() >= OUT_BUF_SIZE) {
ofs[gid] << buf;
buf.clear();
}
}
};

while (true) {
NextRecord nr = reader.getNext();
if (!nr.valid) {
if (haveCurr) flushLine();
break;
}
if (!haveCurr) {
haveCurr = true;
currKmer = nr.kmer;
std::fill(counts.begin(), counts.end(), 0);
}
if (nr.kmer == currKmer) {
counts[nr.genomeID] += nr.count;
} else {
flushLine();
currKmer = nr.kmer;
std::fill(counts.begin(), counts.end(), 0);
counts[nr.genomeID] += nr.count;
}
}

// 5) final flush & close
for (uint8_t gid = 0; gid < nIDs; ++gid) {
if (!buffers[gid].empty()) {
ofs[gid] << buffers[gid];
}
ofs[gid].close();
}

std::cout << "Created per-genome small-k output in " << outputDirectory << "\n";
return 0;
}


int main(int argc, char* argv[]) {
    // Default parameter values
    std::string outputMode = "single";  // Default mode is single 
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

    // read 'k' out of the final.bin header to detect small-k format
    int fileK = 0;
    std::ifstream ifs(finalBin, std::ios::binary);
    char magic[4];
    ifs.read(magic, 4);
    if (std::strncmp(magic, "BINF", 4) != 0) {
        std::cerr << "Bad magic header\n";
        return 1;
    }
    if (!ifs.read(reinterpret_cast<char*>(&fileK), sizeof(fileK))) {
        std::cerr << "Error reading k\n";
        return 1;
    }
    

    // Call the appropriate routine based on the output mode
    if (outputMode == "multiple") {
        if (fileK <= 10) {
            return doMultipleFileSmallKmer(finalBin, outputDirectory);
        } else {
            // existing threaded version for k > 10
            return doMultipleFilesOutput(finalBin, numThreads, outputDirectory);
        }
        
        }   else if (outputMode == "single") {
                    // if k <= 10, use buffered single-threaded dump
                    if (fileK <= 10) {
                        std::string outName = outputFile.empty()
                            ? "final_sorted_" + std::to_string(fileK) + "_dump.txt"
                            : outputFile;
                        return doSingleFileSmallKmer(finalBin, outName);
                    }
                   
                    return doSingleFileOutputWithMmap(finalBin, outputFile, numThreads); //large Ks
                }
}

