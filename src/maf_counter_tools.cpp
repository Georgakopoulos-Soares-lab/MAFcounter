#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <atomic>
#include <thread>
#include <mutex>
#include <cmath>
#include <sstream>
#include <array>
#include <filesystem>
#include <map>
#include <set>
#include <stack>
#include <functional>
#include <cstdint>
#include <cstring>
#include <regex>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
//
// For 128-bit support (boost)
//
#include <boost/multiprecision/cpp_int.hpp>
namespace bmp = boost::multiprecision;
using uint128_t = bmp::uint128_t;

// -------------------------------------------------------------------------------------------
// Some global constants
// -------------------------------------------------------------------------------------------
static const int KMER_MAX_LENGTH = 64;   // Max k
static const int P1_BITS         = 10;   // For big-K partitioning
static const int P2_BITS         = 10;   // For big-K partitioning
static const char BASE_LUT[4]    = {'A','C','G','T'};

// -------------------------------------------------------------------------------------------
// We'll store genome names globally after reading final.bin
// -------------------------------------------------------------------------------------------
static std::vector<std::string> g_genomeNames;
static int      g_k       = 0;    // k read from final.bin
static uint8_t  g_nIDs    = 0;    // number of genome IDs

// -------------------------------------------------------------------------------------------
// Decoding "two-bits-per-base" from a 128-bit integer
// -------------------------------------------------------------------------------------------
static bool g_decodeByteInited = false;
static std::array<std::string, 256> g_decodeByte;

static void initDecodeByte() {
    if (g_decodeByteInited) return;
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
    g_decodeByteInited = true;
}

struct KmerResult {
    std::string kmer;
    uint32_t    count; 
};

// Decode a k-mer from the lower 2*K bits of `combined`.
static KmerResult decodeKmer128(uint128_t combined, int K)
{
    initDecodeByte();

    // Lower 2*K bits for the k-mer
    uint128_t mask = ((uint128_t(1) << (2*K)) - 1);
    uint128_t kmerVal = combined & mask;

    // If the top bits hold a count in your format, you could decode that, but
    // we typically store the count separately in NextRecord, so we do:
    uint32_t cnt = static_cast<uint32_t>((combined >> (2*K)).convert_to<uint64_t>());

    int totalBits = 2*K;
    int leftoverBits = totalBits % 8;
    int fullBytes    = totalBits / 8;

    std::string result;
    result.reserve(K);

    // leftover bits first
    if (leftoverBits > 0) {
        int shift = totalBits - leftoverBits;
        uint64_t partialVal = ((kmerVal >> shift) & ((uint128_t(1) << leftoverBits) - 1)).convert_to<uint64_t>();
        int leftoverBases = leftoverBits / 2;
        for (int i = 0; i < leftoverBases; i++) {
            int baseShift = 2*(leftoverBases-1 - i);
            int baseCode  = (partialVal >> baseShift) & 3;
            result.push_back(BASE_LUT[baseCode]);
        }
    }

    // decode full bytes
    for (int i = fullBytes-1; i >= 0; i--) {
        uint64_t shift   = i*8;
        uint64_t byteVal = ((kmerVal >> shift) & 0xFF).convert_to<uint64_t>();
        result += g_decodeByte[(uint8_t)byteVal];
    }
    if ((int)result.size() > K) {
        result.resize(K);
    }
    return { result, cnt };
}

// -------------------------------------------------------------------------------------------
// Reassemble big-K bits from p1, p2, remainder
// -------------------------------------------------------------------------------------------
static inline uint128_t reassembleKmer(uint16_t p1, uint16_t p2, uint128_t remainder, int k)
{
    // total bits = 2*k
    // p1 => top P1_BITS, p2 => next P2_BITS, remainder => the rest
    int remainder_bits = 2*k - P1_BITS - P2_BITS;
    uint128_t fullVal = p1;
    fullVal <<= (P2_BITS + remainder_bits);
    fullVal |= ((uint128_t)p2 << remainder_bits);
    fullVal |= remainder;
    return fullVal;
}

// Extract the top P1_BITS from the 2*K bits in `fullKmer`
static inline uint16_t getP1Range(uint128_t fullKmer, int k)
{
    int totalBits = 2*k;
    int shift = totalBits - P1_BITS;
    if(shift < 0) return 0;
    uint128_t mask = ((uint128_t(1) << P1_BITS) - 1);
    uint128_t p1val = (fullKmer >> shift) & mask;
    return static_cast<uint16_t>(p1val.convert_to<uint64_t>());
}

// -------------------------------------------------------------------------------------------
// final.bin "big-K block" FileHeader
// -------------------------------------------------------------------------------------------
#pragma pack(push,1)
struct FileHeader {
    char     magic[4];        // "KMHD"
    uint32_t k;
    uint16_t p1_bits;
    uint16_t p2_bits;
    uint8_t  genome_id_bits;
    uint8_t  suffix_mode;     // 0->32bit suffix,1->64bit suffix,3->128bit suffix
    uint8_t  reserved[2];
};
#pragma pack(pop)

// A structure for reading next record from final.bin in ascending k-mer order
struct NextRecord {
    bool      valid;   
    uint128_t kmer;    
    uint8_t   genomeID;
    uint32_t  count;   
};

// -------------------------------------------------------------------------------------------
// A class that streams from final.bin in ascending k-mer order, restricted to p1 range
// -------------------------------------------------------------------------------------------
class FinalBinReader {
public:
    FinalBinReader(const std::string &binPath,
                   const std::string &metaPath,
                   int k,
                   uint16_t p1Start,
                   uint16_t p1End)
        : m_k(k), m_p1Start(p1Start), m_p1End(p1End)
    {
        m_ifs.open(binPath, std::ios::binary);
        if(!m_ifs) {
            m_eof = true;
            std::cerr<<"Cannot open "<<binPath<<"\n";
            return;
        }
        // read "BINF"
        {
            char mg[4];
            if(!m_ifs.read(mg,4)){m_eof=true;return;}
            if(std::strncmp(mg,"BINF",4)!=0){
                std::cerr<<"Not a valid BINF\n";
                m_eof=true;return;
            }
        }
        // read stored k (not necessarily cross-check with user k)
        {
            int fileK=0;
            if(!m_ifs.read((char*)&fileK,sizeof(fileK))){m_eof=true;return;}
        }
        // read nIDs + skip genome names
        {
            uint8_t nIDs=0;
            if(!m_ifs.read((char*)&nIDs,sizeof(nIDs))){m_eof=true;return;}
            for(uint8_t i=0; i<nIDs; i++){
                uint16_t len=0;
                if(!m_ifs.read((char*)&len,sizeof(len))){m_eof=true;return;}
                m_ifs.ignore(len);
            }
        }
        // read "ENDG"
        {
            char sentinel[4];
            if(!m_ifs.read(sentinel,4)){m_eof=true;return;}
            if(std::strncmp(sentinel,"ENDG",4)!=0){
                std::cerr<<"Missing ENDG sentinel\n";
                m_eof=true;return;
            }
        }

        // small-K vs big-K
        m_smallKmode = (m_k <= 10);

        if(m_smallKmode) {
            // read mapSize
            if(!m_ifs.read((char*)&m_mapSize,sizeof(m_mapSize))){m_eof=true;return;}
            m_mapReadCount=0;
            m_idx=256;
        } else {
            // for big-K => we need metadata
            std::ifstream meta(metaPath);
            if(!meta){
                std::cerr<<"Cannot open metadata file: "<<metaPath<<"\n";
                m_eof=true;return;
            }
            std::vector<std::pair<uint16_t,uint64_t>> metaEntries;
            std::string line;
            while(std::getline(meta,line)){
                if(line.empty()) continue;
                std::istringstream iss(line);
                uint16_t p1val;
                uint64_t off;
                if(!(iss >> p1val >> off)) continue;
                metaEntries.push_back({p1val, off});
            }
            meta.close();
            if(metaEntries.empty()){
                std::cerr<<"Metadata is empty?\n";
                m_eof=true;return;
            }
            std::sort(metaEntries.begin(), metaEntries.end(),
                      [](auto &a, auto &b){return a.first < b.first;});
            auto it = std::lower_bound(metaEntries.begin(), metaEntries.end(), p1Start,
                [](const std::pair<uint16_t,uint64_t> &e, uint16_t val){
                    return e.first < val;
                }
            );
            if(it==metaEntries.end()){
                m_eof=true;return;
            }
            m_ifs.seekg(it->second);
            if(!m_ifs){
                m_eof=true;return;
            }
            m_inBinBlock=false;
        }
    }

    bool good() const { return !m_eof; }

    NextRecord getNext()
    {
        if(m_eof) return {false,0,0,0};
        if(m_smallKmode) return getNextSmall();
        else return getNextBig();
    }

private:
    NextRecord getNextSmall()
    {
        while(true){
            if(m_eof) return {false,0,0,0};
            if(m_mapReadCount>=m_mapSize){
                m_eof=true;
                return {false,0,0,0};
            }
            if(m_idx>=256){
                // read key + 256 counts
                uint32_t key32=0;
                if(!m_ifs.read((char*)&key32,sizeof(key32))){
                    m_eof=true;return {false,0,0,0};
                }
                if(!m_ifs.read((char*)m_counts.data(),256*sizeof(uint32_t))){
                    m_eof=true;return {false,0,0,0};
                }
                m_mapReadCount++;
                m_idx=0;
                m_currKmer = key32;  // fits into bottom bits of uint128_t
            }
            while(m_idx<256 && m_counts[m_idx]==0){
                m_idx++;
            }
            if(m_idx>=256) continue; // read next entry
            NextRecord nr;
            nr.valid    = true;
            nr.kmer     = m_currKmer;
            nr.genomeID = (uint8_t)m_idx;
            nr.count    = m_counts[m_idx];
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
            if(m_eof) return {false,0,0,0};
            if(!m_inBinBlock){
                // read "BINc"
                char binC[4];
                if(!m_ifs.read(binC,4)){m_eof=true;return{false,0,0,0};}
                if(std::strncmp(binC,"BINc",4)!=0){
                    m_eof=true;return{false,0,0,0};
                }
                if(!m_ifs.read((char*)&m_blockP1,sizeof(m_blockP1))){
                    m_eof=true;return{false,0,0,0};
                }
                if(m_blockP1>m_p1End){
                    m_eof=true;return{false,0,0,0};
                }
                // read "KMHD"
                char kmhd[4];
                if(!m_ifs.read(kmhd,4)){m_eof=true;return{false,0,0,0};}
                if(std::strncmp(kmhd,"KMHD",4)!=0){
                    m_eof=true;return{false,0,0,0};
                }
                FileHeader hdr;
                std::memcpy(hdr.magic, kmhd,4);
                if(!m_ifs.read((char*)&hdr.k,sizeof(hdr.k))){
                    m_eof=true;return{false,0,0,0};
                }
                if(!m_ifs.read((char*)&hdr.p1_bits,sizeof(hdr.p1_bits))){
                    m_eof=true;return{false,0,0,0};
                }
                if(!m_ifs.read((char*)&hdr.p2_bits,sizeof(hdr.p2_bits))){
                    m_eof=true;return{false,0,0,0};
                }
                if(!m_ifs.read((char*)&hdr.genome_id_bits,sizeof(hdr.genome_id_bits))){
                    m_eof=true;return{false,0,0,0};
                }
                if(!m_ifs.read((char*)&hdr.suffix_mode,sizeof(hdr.suffix_mode))){
                    m_eof=true;return{false,0,0,0};
                }
                if(!m_ifs.read((char*)&hdr.reserved[0],2)){
                    m_eof=true;return{false,0,0,0};
                }
                // read "GRPS"
                char s1[4];
                if(!m_ifs.read(s1,4)){m_eof=true;return{false,0,0,0};}
                uint32_t groupCount=0;
                if(!m_ifs.read((char*)&groupCount,sizeof(groupCount))){
                    m_eof=true;return{false,0,0,0};
                }
                m_groups.clear();
                m_groups.resize(groupCount);
                for(uint32_t i=0;i<groupCount;i++){
                    if(!m_ifs.read((char*)&m_groups[i].p2,sizeof(m_groups[i].p2))){
                        m_eof=true;return{false,0,0,0};
                    }
                    if(!m_ifs.read((char*)&m_groups[i].sz,sizeof(m_groups[i].sz))){
                        m_eof=true;return{false,0,0,0};
                    }
                }
                // read "ENDS"
                char s2[4];
                if(!m_ifs.read(s2,4)){m_eof=true;return{false,0,0,0};}
                m_suffixMode = hdr.suffix_mode;
                m_inBinBlock=true;
                m_groupIndex=0;
            }
            // skip done groups
            while(m_groupIndex < m_groups.size() && m_groups[m_groupIndex].done()){
                m_groupIndex++;
            }
            if(m_groupIndex>=m_groups.size()){
                // read "BEND"
                char bend[4];
                if(!m_ifs.read(bend,4)){m_eof=true;return{false,0,0,0};}
                m_inBinBlock=false;
                continue;
            }
            auto &cg = m_groups[m_groupIndex];
            if(!cg.loaded){
                if(m_suffixMode==0){
                    cg.data32.resize(cg.sz);
                    if(!m_ifs.read((char*)cg.data32.data(), cg.sz*sizeof(uint32_t))){
                        m_eof=true;return{false,0,0,0};
                    }
                } else if(m_suffixMode==1){
                    cg.data64.resize(cg.sz);
                    if(!m_ifs.read((char*)cg.data64.data(), cg.sz*sizeof(uint64_t))){
                        m_eof=true;return{false,0,0,0};
                    }
                } else if(m_suffixMode==3){
                    cg.data128.resize(cg.sz);
                    if(!m_ifs.read((char*)cg.data128.data(), cg.sz*2*sizeof(uint64_t))){
                        m_eof=true;return{false,0,0,0};
                    }
                } else {
                    std::cerr<<"Unknown suffix_mode "<<(int)m_suffixMode<<"\n";
                    m_eof=true;return{false,0,0,0};
                }
                cg.loaded=true;
                cg.idx=0;
            }
            // read next from group
            if(m_suffixMode==0){
                if(cg.idx>=cg.data32.size()){
                    cg.idx=cg.sz; 
                    continue;
                }
                uint32_t val = cg.data32[cg.idx++];
                uint8_t gID  = (uint8_t)(val & 0xFF);
                uint64_t remainder = (val >> 8);
                uint128_t fullK = reassembleKmer(m_blockP1,cg.p2,remainder,m_k);
                uint16_t p1 = getP1Range(fullK,m_k);
                if(p1>=m_p1Start && p1<=m_p1End){
                    return {true,fullK,gID,1};
                }
            } else if(m_suffixMode==1){
                if(cg.idx>=cg.data64.size()){
                    cg.idx=cg.sz;
                    continue;
                }
                uint64_t val = cg.data64[cg.idx++];
                uint8_t gID  = (uint8_t)(val & 0xFF);
                uint64_t remainder = (val >> 8);
                uint128_t fullK = reassembleKmer(m_blockP1,cg.p2,remainder,m_k);
                uint16_t p1 = getP1Range(fullK,m_k);
                if(p1>=m_p1Start && p1<=m_p1End){
                    return {true,fullK,gID,1};
                }
            } else {
                // suffixMode==3 => 128 bits
                if(cg.idx>=cg.data128.size()){
                    cg.idx=cg.sz;
                    continue;
                }
                auto &arr = cg.data128[cg.idx++];
                uint8_t gID = (uint8_t)(arr[0] & 0xFF);
                uint128_t remainder = (((uint128_t)arr[1]) << 56)
                                      | ((arr[0] >> 8) & 0x00FFFFFFFFFFFFFFULL);
                uint128_t fullK = reassembleKmer(m_blockP1,cg.p2,remainder,m_k);
                uint16_t p1 = getP1Range(fullK,m_k);
                if(p1>=m_p1Start && p1<=m_p1End){
                    return {true,fullK,gID,1};
                }
            }
        }
    }

private:
    std::ifstream m_ifs;
    bool m_eof=false;
    bool m_smallKmode=false;
    int  m_k=0;

    // for small-K
    uint64_t m_mapSize=0;
    uint64_t m_mapReadCount=0;
    std::array<uint32_t,256> m_counts{};
    uint128_t m_currKmer=0;
    int m_idx=256;

    // for big-K
    bool m_inBinBlock=false;
    uint16_t m_blockP1=0;
    uint8_t  m_suffixMode=0;

    struct GroupData {
        uint16_t p2=0;
        uint32_t sz=0;
        bool loaded=false;
        std::vector<uint32_t> data32;
        std::vector<uint64_t> data64;
        std::vector<std::array<uint64_t,2>> data128;
        size_t idx=0;
        bool done() const {return idx>=sz;}
    };
    std::vector<GroupData> m_groups;
    size_t m_groupIndex=0;

    uint16_t m_p1Start=0;
    uint16_t m_p1End=0;
};

// -------------------------------------------------------------------------------------------
// Utility: compute standard deviation for non-zero counts
// -------------------------------------------------------------------------------------------
static double computeStdDev(const std::vector<uint32_t> &counts)
{
    std::vector<double> nonZero;
    nonZero.reserve(counts.size());
    for(uint32_t c : counts){
        if(c>0) nonZero.push_back((double)c);
    }
    if(nonZero.size()<2){
        return 0.0;
    }
    double sum=0;
    for(auto d: nonZero) sum+=d;
    double mean = sum / nonZero.size();
    double sq=0;
    for(auto d: nonZero){
        double diff=(d-mean);
        sq += diff*diff;
    }
    double variance = sq / nonZero.size();
    return std::sqrt(variance);
}

// -------------------------------------------------------------------------------------------
// Data structure for top-std mode
// -------------------------------------------------------------------------------------------
struct TopStdRecord {
    double stddev;
    uint128_t kmer;
    std::vector<uint32_t> counts;
};

static std::vector<TopStdRecord> g_threadTopResults; // for doComputeTopStds
struct ThreadRange {
    uint16_t p1Start;
    uint16_t p1End;
};

static int g_topCount=0;       
static int g_numThreads=1;     

// Worker thread for --std
static void workerThread_std(int threadID,
                             ThreadRange range,
                             const std::string &binPath,
                             const std::string &metaPath)
{
    FinalBinReader reader(binPath, metaPath, g_k, range.p1Start, range.p1End);
    if(!reader.good()){
        // fill with zero
        for(int i=0; i<g_topCount;i++){
            g_threadTopResults[threadID*g_topCount + i] = {0.0,0,{}};
        }
        return;
    }
    auto cmpMinHeap = [](const TopStdRecord &a, const TopStdRecord &b){
        return a.stddev > b.stddev; 
    };
    std::vector<TopStdRecord> heap;
    heap.reserve(g_topCount+1);

    bool haveCurr=false;
    uint128_t currKmer=0;
    std::vector<uint32_t> counts(g_nIDs,0);

    auto flush=[&](){
        double sd = computeStdDev(counts);
        if(sd>0.0){
            TopStdRecord rec;
            rec.stddev=sd;
            rec.kmer=currKmer;
            rec.counts=counts;
            if((int)heap.size()<g_topCount){
                heap.push_back(std::move(rec));
                std::push_heap(heap.begin(), heap.end(), cmpMinHeap);
            } else {
                if(rec.stddev > heap.front().stddev){
                    std::pop_heap(heap.begin(), heap.end(), cmpMinHeap);
                    heap.pop_back();
                    heap.push_back(std::move(rec));
                    std::push_heap(heap.begin(), heap.end(), cmpMinHeap);
                }
            }
        }
    };

    while(true){
        NextRecord nr= reader.getNext();
        if(!nr.valid){
            if(haveCurr) flush();
            break;
        }
        if(!haveCurr){
            haveCurr=true;
            currKmer= nr.kmer;
            std::fill(counts.begin(), counts.end(), 0);
            counts[nr.genomeID]+= nr.count;
        } else {
            if(nr.kmer==currKmer){
                counts[nr.genomeID]+= nr.count;
            } else {
                flush();
                currKmer= nr.kmer;
                std::fill(counts.begin(), counts.end(), 0);
                counts[nr.genomeID]+= nr.count;
            }
        }
    }

    std::sort_heap(heap.begin(), heap.end(), cmpMinHeap);
    int base=threadID*g_topCount;
    for(int i=0; i<(int)heap.size(); i++){
        g_threadTopResults[base + i] = std::move(heap[i]);
    }
    for(int i=(int)heap.size(); i<g_topCount;i++){
        g_threadTopResults[base+i]={0.0,0,{}};
    }
}

// doComputeTopStds
static int doComputeTopStds(const std::string &binPath,
                            const std::string &metaPath,
                            int topCount,
                            int numThreads)
{
    // read final.bin header
    {
        std::ifstream ifs(binPath,std::ios::binary);
        if(!ifs){
            std::cerr<<"Cannot open "<<binPath<<"\n";
            return 1;
        }
        {
            char mg[4];
            if(!ifs.read(mg,4)){
                std::cerr<<"Cannot read magic\n";return 1;
            }
            if(std::strncmp(mg,"BINF",4)!=0){
                std::cerr<<"Not a valid BINF\n";return 1;
            }
        }
        {
            if(!ifs.read((char*)&g_k,sizeof(g_k))){
                std::cerr<<"Cannot read k\n";return 1;
            }
            if(g_k<1||g_k> KMER_MAX_LENGTH){
                std::cerr<<"Invalid k\n";return 1;
            }
        }
        {
            if(!ifs.read((char*)&g_nIDs,sizeof(g_nIDs))){
                std::cerr<<"Cannot read nIDs\n";return 1;
            }
            if(g_nIDs<1){
                std::cerr<<"No genome IDs?\n";return 1;
            }
        }
        g_genomeNames.resize(g_nIDs);
        for(uint8_t i=0;i<g_nIDs;i++){
            uint16_t len=0;
            if(!ifs.read((char*)&len,sizeof(len))){
                std::cerr<<"Cannot read genome name length\n";return 1;
            }
            g_genomeNames[i].resize(len);
            if(!ifs.read((char*)g_genomeNames[i].data(),len)){
                std::cerr<<"Cannot read genome name\n";return 1;
            }
        }
        {
            char sentinel[4];
            if(!ifs.read(sentinel,4)){
                std::cerr<<"Cannot read ENDG\n";return 1;
            }
            if(std::strncmp(sentinel,"ENDG",4)!=0){
                std::cerr<<"Missing ENDG\n";return 1;
            }
        }
    }

    g_topCount = topCount;
    g_numThreads = numThreads;
    g_threadTopResults.resize(g_topCount*g_numThreads);

    int totalP1 = (1<<P1_BITS);
    int chunkSize = (totalP1 + g_numThreads -1)/g_numThreads;
    std::vector<ThreadRange> ranges(g_numThreads);
    for(int i=0;i<g_numThreads;i++){
        uint16_t st = (uint16_t)(i*chunkSize);
        uint16_t ed = (uint16_t)std::min<int>((i+1)*chunkSize-1, totalP1-1);
        ranges[i]={st,ed};
    }

    std::vector<std::thread> threads;
    for(int i=0;i<g_numThreads;i++){
        threads.emplace_back(workerThread_std, i, ranges[i], binPath, metaPath);
    }
    for(auto &t: threads){
        t.join();
    }
    // gather top
    std::sort(g_threadTopResults.begin(), g_threadTopResults.end(),
              [](auto &a, auto &b){return a.stddev > b.stddev;});
    if((int)g_threadTopResults.size()>g_topCount){
        g_threadTopResults.resize(g_topCount);
    }

    std::string outName = "top_"+std::to_string(g_topCount)+"_std.txt";
    std::ofstream ofs(outName);
    if(!ofs){
        std::cerr<<"Cannot create "<<outName<<"\n";
        return 1;
    }
    ofs << "# kmer, stddev, counts\n";
    for(auto &rec: g_threadTopResults){
        if(rec.stddev<=0.0) break;
        KmerResult kr = decodeKmer128(rec.kmer, g_k);
        ofs<< kr.kmer << ", " << rec.stddev << ", ";
        bool first=true;
        for(uint8_t gid=0; gid<g_nIDs; gid++){
            uint32_t c= rec.counts[gid];
            if(c>0){
                if(!first) ofs<<"; ";
                first=false;
                if(gid<g_genomeNames.size()){
                    ofs << g_genomeNames[gid] <<": "<< c;
                } else {
                    ofs << "UnknownID#"<<(int)gid <<": "<< c;
                }
            }
        }
        ofs<<"\n";
    }
    ofs.close();
    std::cout<<"Wrote "<<outName<<" with up to "<<g_topCount<<" highest-stddev k-mers.\n";
    return 0;
}

// -------------------------------------------------------------------------------------------
// Expression parsing for --expr
// -------------------------------------------------------------------------------------------
enum class TokenType {
    VAR,
    OP,
    AND,
    OR,
    LPAREN,
    RPAREN,
    NUMBER,
    END
};

struct Token {
    TokenType type;
    std::string text;
};

static std::vector<Token> tokenizeExpr(const std::string &expr)
{
    std::vector<Token> tokens;
    size_t i=0, n=expr.size();
    auto skipSpaces=[&](){
        while(i<n && std::isspace((unsigned char)expr[i])) i++;
    };
    skipSpaces();
    while(i<n){
        if(std::isspace((unsigned char)expr[i])){
            skipSpaces();
            continue;
        }
        if(expr[i]=='('){
            tokens.push_back({TokenType::LPAREN,"("});
            i++;
            continue;
        } else if(expr[i]==')'){
            tokens.push_back({TokenType::RPAREN,")"});
            i++;
            continue;
        }
        if(i+1<n){
            if(expr[i]=='&' && expr[i+1]=='&'){
                tokens.push_back({TokenType::AND,"&&"});
                i+=2;
                continue;
            }
            if(expr[i]=='|' && expr[i+1]=='|'){
                tokens.push_back({TokenType::OR,"||"});
                i+=2;
                continue;
            }
        }
        if(expr[i]=='>' || expr[i]=='<' || expr[i]=='='){
            char c=expr[i];
            if(i+1<n && expr[i+1]=='='){
                std::string t; t.push_back(c); t.push_back('=');
                tokens.push_back({TokenType::OP,t});
                i+=2;
            } else {
                std::string t(1,c);
                tokens.push_back({TokenType::OP,t});
                i++;
            }
            continue;
        }
        // digit => parse number
        if(std::isdigit((unsigned char)expr[i])){
            size_t start=i;
            while(i<n && std::isdigit((unsigned char)expr[i])) i++;
            tokens.push_back({TokenType::NUMBER, expr.substr(start,i-start)});
            continue;
        }
        // parse variable
        {
            size_t start=i;
            while(i<n && !std::isspace((unsigned char)expr[i]) &&
                  expr[i]!=')' && expr[i]!='(' &&
                  !(i+1<n && ((expr[i]=='&' && expr[i+1]=='&')|| (expr[i]=='|' && expr[i+1]=='|'))) &&
                  expr[i]!='<' && expr[i]!='>' && expr[i]!='='
                 ){
                i++;
            }
            std::string var = expr.substr(start, i-start);
            tokens.push_back({TokenType::VAR,var});
        }
    }
    tokens.push_back({TokenType::END,""});
    return tokens;
}

struct Condition {
    int genomeIndex;
    std::string op;
    int value;
};

static bool evalCondition(const Condition &cond, const std::vector<uint32_t> &counts)
{
    uint32_t c= counts[cond.genomeIndex];
    if(cond.op==">")  return (c> (uint32_t)cond.value);
    if(cond.op==">=") return (c>=(uint32_t)cond.value);
    if(cond.op=="<")  return (c< (uint32_t)cond.value);
    if(cond.op=="<=") return (c<=(uint32_t)cond.value);
    if(cond.op=="=")  return (c== (uint32_t)cond.value);
    return false;
}

static std::function<bool(const std::vector<uint32_t>&)>
buildExprEvaluator(const std::string &exprStr,
                   const std::unordered_map<std::string,int> &varIndex)
{
    auto tokens = tokenizeExpr(exprStr);

    struct Item {
        bool isCondition=false;
        Condition cond;
        std::string op; 
    };

    auto parseCondition=[&](size_t &idx)->Condition {
        if(idx+2>=tokens.size()){
            throw std::runtime_error("Incomplete condition near expression end");
        }
        if(tokens[idx].type!=TokenType::VAR ||
           tokens[idx+1].type!=TokenType::OP ||
           tokens[idx+2].type!=TokenType::NUMBER){
            throw std::runtime_error("Expected: VAR OP NUMBER");
        }
        std::string varName = tokens[idx].text;
        std::string op      = tokens[idx+1].text;
        int val             = std::stoi(tokens[idx+2].text);

        auto it= varIndex.find(varName);
        if(it==varIndex.end()){
            throw std::runtime_error("Variable '"+varName+"' not recognized among genome names");
        }
        Condition c;
        c.genomeIndex= it->second;
        c.op=op;
        c.value= val;
        idx+=3;
        return c;
    };

    // Shunting-yard
    std::vector<Item> output;
    std::stack<std::string> opStack;

    auto precedenceOf=[&](const std::string &op)->int{
        if(op=="AND") return 2;
        if(op=="OR")  return 1;
        return 0;
    };

    size_t i=0;
    while(true){
        if(tokens[i].type==TokenType::END) break;
        if(tokens[i].type==TokenType::VAR){
            Condition c = parseCondition(i);
            output.push_back({true,c,""});
        }
        else if(tokens[i].type==TokenType::LPAREN){
            opStack.push("(");
            i++;
        }
        else if(tokens[i].type==TokenType::RPAREN){
            while(!opStack.empty() && opStack.top()!="("){
                std::string topOp= opStack.top();
                opStack.pop();
                output.push_back({false,{},topOp});
            }
            if(opStack.empty()){
                throw std::runtime_error("Mismatched parentheses");
            }
            opStack.pop();
            i++;
        }
        else if(tokens[i].type==TokenType::AND || tokens[i].type==TokenType::OR){
            std::string thisOp = (tokens[i].type==TokenType::AND) ? "AND":"OR";
            i++;
            while(!opStack.empty() && opStack.top()!="(" &&
                  precedenceOf(opStack.top())>=precedenceOf(thisOp)){
                output.push_back({false,{},opStack.top()});
                opStack.pop();
            }
            opStack.push(thisOp);
        }
        else {
            throw std::runtime_error("Syntax error near: "+tokens[i].text);
        }
    }
    while(!opStack.empty()){
        if(opStack.top()=="("){
            throw std::runtime_error("Mismatched parentheses in expression");
        }
        output.push_back({false,{},opStack.top()});
        opStack.pop();
    }

    // build a lambda for evaluation
    return [output](const std::vector<uint32_t> &counts)->bool {
        std::stack<bool> st;
        for(auto &item: output){
            if(item.isCondition){
                bool b= evalCondition(item.cond, counts);
                st.push(b);
            } else {
                // AND/OR
                if(st.size()<2) return false;
                bool b2=st.top(); st.pop();
                bool b1=st.top(); st.pop();
                if(item.op=="AND") st.push(b1 && b2);
                else st.push(b1 || b2);
            }
        }
        if(st.size()!=1) return false;
        return st.top();
    };
}

// -------------------------------------------------------------------------------------------
// --expr => 2 pass approach => "filtered_expr.txt" with matching k-mers
// -------------------------------------------------------------------------------------------


static std::array<std::atomic<size_t>, 1024> g_p1SizesAtomic; 
static std::array<size_t, 1025> g_p1Offsets;

static int digitCount(uint32_t num){
    if(num==0) return 1;
    int d=0;
    while(num>0){
        num/=10;
        d++;
    }
    return d;
}

static size_t computeLineSize(const std::string &kmerStr,
                              const std::unordered_map<uint8_t,uint32_t> &gcounts)
{
    // Format: "kmer: name: count, name2: count2\n"
    // We'll assume the genome name is stored in g_genomeNames
    size_t lineSize=0;
    // kmer + ": "
    lineSize += kmerStr.size()+2;
    bool first=true;
    for(auto &kv: gcounts){
        if(kv.second==0) continue;
        if(!first) lineSize+=2; // ", "
        first=false;
        uint8_t gID= kv.first;
        if(gID<g_genomeNames.size()){
            lineSize+= g_genomeNames[gID].size(); 
        } else {
            static const char pfx[]="UnknownID#";
            lineSize += 9 + digitCount(gID);
        }
        lineSize+=2; // ": "
        lineSize+= digitCount(kv.second);
    }
    lineSize+=1; // newline
    return lineSize;
}

// pass1
static void p1SizeComputerWorker_expr(int threadID,
                                      ThreadRange range,
                                      const std::string &binPath,
                                      const std::string &metaPath,
                                      const std::function<bool(const std::vector<uint32_t>&)> &exprFn)
{
    FinalBinReader reader(binPath, metaPath, g_k, range.p1Start, range.p1End);
    if(!reader.good()) return;
    bool haveCurr=false;
    uint128_t currKmer=0;
    std::vector<uint32_t> counts(g_nIDs,0);

    auto flush=[&](){
        if(!exprFn(counts)) return;
        std::unordered_map<uint8_t,uint32_t> gmap;
        for(uint8_t i=0; i<g_nIDs; i++){
            if(counts[i]>0) gmap[i]=counts[i];
        }
        if(gmap.empty()) return;
        KmerResult kr = decodeKmer128(currKmer, g_k);
        uint16_t p1= getP1Range(currKmer, g_k);
        size_t inc= computeLineSize(kr.kmer, gmap);
        g_p1SizesAtomic[p1].fetch_add(inc, std::memory_order_relaxed);
    };

    while(true){
        NextRecord nr=reader.getNext();
        if(!nr.valid){
            if(haveCurr) flush();
            break;
        }
        if(!haveCurr){
            haveCurr=true;
            currKmer=nr.kmer;
            std::fill(counts.begin(),counts.end(),0);
            counts[nr.genomeID]+=nr.count;
        } else {
            if(nr.kmer==currKmer){
                counts[nr.genomeID]+= nr.count;
            } else {
                flush();
                currKmer= nr.kmer;
                std::fill(counts.begin(),counts.end(),0);
                counts[nr.genomeID]+= nr.count;
            }
        }
    }
}

static void buildLine(char* dest,
                      const std::string &kmerStr,
                      const std::unordered_map<uint8_t,uint32_t> &gmap,
                      size_t &written)
{
    memcpy(dest+written, kmerStr.data(), kmerStr.size());
    written+= kmerStr.size();
    dest[written++]=':';
    dest[written++]=' ';
    bool first=true;
    for(auto &kv: gmap){
        if(!first){
            dest[written++]=',';
            dest[written++]=' ';
        }
        first=false;
        uint8_t gID= kv.first;
        if(gID< g_genomeNames.size()){
            const std::string &nm= g_genomeNames[gID];
            memcpy(dest+written, nm.data(), nm.size());
            written+= nm.size();
        } else {
            static const char pfx[]="UnknownID#";
            memcpy(dest+written, pfx, 9);
            written+=9;
            char buf[32];
            int len= std::snprintf(buf,32,"%u",(unsigned)gID);
            memcpy(dest+written, buf, len);
            written+=len;
        }
        dest[written++]=':';
        dest[written++]=' ';
        char cbuf[32];
        int len= std::snprintf(cbuf,32,"%u",kv.second);
        memcpy(dest+written, cbuf, len);
        written+=len;
    }
    dest[written++]='\n';
}

// pass2
static void p1WriteWorker_expr(int threadID,
                               ThreadRange range,
                               const std::string &binPath,
                               const std::string &metaPath,
                               const std::function<bool(const std::vector<uint32_t>&)> &exprFn,
                               char* mappedPtr)
{
    FinalBinReader reader(binPath, metaPath, g_k, range.p1Start, range.p1End);
    if(!reader.good()) return;

    bool haveCurr=false;
    uint128_t currKmer=0;
    std::vector<uint32_t> counts(g_nIDs,0);

    auto flush=[&](){
        if(!exprFn(counts)) return;
        std::unordered_map<uint8_t,uint32_t> gmap;
        for(uint8_t i=0; i<g_nIDs; i++){
            if(counts[i]>0) gmap[i]=counts[i];
        }
        if(gmap.empty()) return;
        KmerResult kr= decodeKmer128(currKmer,g_k);
        uint16_t p1= getP1Range(currKmer,g_k);
        size_t offset= g_p1Offsets[p1];
        buildLine(mappedPtr, kr.kmer, gmap, offset);
        g_p1Offsets[p1]= offset;
    };

    while(true){
        NextRecord nr= reader.getNext();
        if(!nr.valid){
            if(haveCurr) flush();
            break;
        }
        if(!haveCurr){
            haveCurr=true;
            currKmer= nr.kmer;
            std::fill(counts.begin(), counts.end(),0);
            counts[nr.genomeID]+= nr.count;
        } else {
            if(nr.kmer==currKmer){
                counts[nr.genomeID]+= nr.count;
            } else {
                flush();
                currKmer= nr.kmer;
                std::fill(counts.begin(),counts.end(),0);
                counts[nr.genomeID]+= nr.count;
            }
        }
    }
}

static int doComputeExpr(const std::string &binPath,
                         const std::string &metaPath,
                         int numThreads,
                         const std::string &exprStr)
{
    // read final.bin header
    {
        std::ifstream ifs(binPath,std::ios::binary);
        if(!ifs){
            std::cerr<<"Cannot open "<<binPath<<"\n";
            return 1;
        }
        char mg[4];
        if(!ifs.read(mg,4)){
            std::cerr<<"Cannot read magic\n";return 1;
        }
        if(std::strncmp(mg,"BINF",4)!=0){
            std::cerr<<"Not a valid BINF\n";return 1;
        }
        if(!ifs.read((char*)&g_k,sizeof(g_k))){
            std::cerr<<"Cannot read k\n";return 1;
        }
        if(g_k<1 || g_k> KMER_MAX_LENGTH){
            std::cerr<<"Invalid k\n";return 1;
        }
        if(!ifs.read((char*)&g_nIDs,sizeof(g_nIDs))){
            std::cerr<<"Cannot read nIDs\n";return 1;
        }
        if(g_nIDs<1){
            std::cerr<<"No genome IDs?\n";return 1;
        }
        g_genomeNames.resize(g_nIDs);
        for(uint8_t i=0;i<g_nIDs;i++){
            uint16_t len=0;
            if(!ifs.read((char*)&len,sizeof(len))){
                std::cerr<<"Cannot read genome name length\n";return 1;
            }
            g_genomeNames[i].resize(len);
            if(!ifs.read((char*)g_genomeNames[i].data(),len)){
                std::cerr<<"Cannot read genome name\n";return 1;
            }
        }
        {
            char sentinel[4];
            if(!ifs.read(sentinel,4)){
                std::cerr<<"Cannot read ENDG\n";return 1;
            }
            if(std::strncmp(sentinel,"ENDG",4)!=0){
                std::cerr<<"Missing ENDG\n";return 1;
            }
        }
    }

    // build varIndex
    std::unordered_map<std::string,int> varIndex;
    for(int i=0;i<(int)g_genomeNames.size(); i++){
        varIndex[g_genomeNames[i]]= i;
    }
    std::function<bool(const std::vector<uint32_t>&)> exprFn;
    try {
        exprFn= buildExprEvaluator(exprStr, varIndex);
    } catch(std::exception &ex){
        std::cerr<<"Expression parse error: "<<ex.what()<<"\n";
        return 1;
    }

    for(auto &x: g_p1SizesAtomic){
        x.store(0, std::memory_order_relaxed);
    }
    g_numThreads= std::max(numThreads,1);

    int totalP1=(1<<P1_BITS);
    int chunkSize=(totalP1+g_numThreads-1)/g_numThreads;
    std::vector<ThreadRange> thrRanges(g_numThreads);
    for(int i=0;i<g_numThreads;i++){
        uint16_t st=(uint16_t)(i*chunkSize);
        uint16_t ed=(uint16_t)std::min<int>((i+1)*chunkSize-1,totalP1-1);
        thrRanges[i]={st,ed};
    }

    // pass1
    {
        std::vector<std::thread> threads;
        for(int t=0;t<g_numThreads;t++){
            threads.emplace_back(p1SizeComputerWorker_expr,
                                 t, thrRanges[t], binPath, metaPath,
                                 std::cref(exprFn));
        }
        for(auto &th: threads){
            th.join();
        }
    }

    // prefix sum
    size_t totalBytes=0;
    for(int i=0;i<1024;i++){
        size_t sz= g_p1SizesAtomic[i].load(std::memory_order_relaxed);
        g_p1Offsets[i]= totalBytes;
        totalBytes += sz;
    }
    g_p1Offsets[1024]= totalBytes;

    std::string outName="filtered_expr.txt";
    int fd= ::open(outName.c_str(), O_RDWR|O_CREAT, 0666);
    if(fd<0){
        std::cerr<<"Cannot create "<<outName<<"\n";
        return 1;
    }
    if(ftruncate(fd,(off_t)totalBytes)!=0){
        std::cerr<<"ftruncate failed\n";
        close(fd);
        return 1;
    }
    void* mapPtr= mmap(nullptr,totalBytes,PROT_READ|PROT_WRITE, MAP_SHARED, fd, 0);
    if(mapPtr==MAP_FAILED){
        std::cerr<<"mmap failed\n";
        close(fd);
        return 1;
    }
    char* mappedPtr=(char*)mapPtr;

    // pass2
    {
        std::vector<std::thread> threads;
        for(int t=0;t<g_numThreads;t++){
            threads.emplace_back(p1WriteWorker_expr,
                                 t, thrRanges[t],
                                 binPath, metaPath,
                                 std::cref(exprFn),
                                 mappedPtr);
        }
        for(auto &th: threads){
            th.join();
        }
    }

    msync(mappedPtr, totalBytes, MS_SYNC);
    munmap(mappedPtr,totalBytes);
    close(fd);
    std::cout<<"Created "<<outName<<" (size "<<totalBytes<<" bytes) with expression-filtered k-mers.\n";
    return 0;
}

// -------------------------------------------------------------------------------------------
// --gstats => compute distribution stats per genome
// -------------------------------------------------------------------------------------------
struct PartialGstats {
    std::vector<uint32_t> minCount;
    std::vector<uint32_t> maxCount;
    std::vector<uint128_t> sum1;
    std::vector<uint128_t> sum2;
    std::vector<uint128_t> sum3;
    std::vector<uint64_t>  nKmers;
    std::vector<std::unordered_map<uint32_t,uint64_t>> freqMap;
};

static std::vector<PartialGstats> g_threadPartials;

static void initPartialGstats(PartialGstats &p, int nIDs)
{
    p.minCount.resize(nIDs,(uint32_t)-1);
    p.maxCount.resize(nIDs,0);
    p.sum1.resize(nIDs,0);
    p.sum2.resize(nIDs,0);
    p.sum3.resize(nIDs,0);
    p.nKmers.resize(nIDs,0);
    p.freqMap.resize(nIDs);
}

static void updatePartial(PartialGstats &p, const std::vector<uint32_t> &counts)
{
    int n = (int)counts.size();
    for(int i=0;i<n;i++){
        uint32_t c= counts[i];
        if(c< p.minCount[i]) p.minCount[i]= c;
        if(c> p.maxCount[i]) p.maxCount[i]= c;
        p.sum1[i]+= c;
        uint128_t c2=(uint128_t)c*c;
        p.sum2[i]+= c2;
        p.sum3[i]+= (c2*c);
        p.nKmers[i]+=1;
        p.freqMap[i][c]+=1;
    }
}

static void workerThread_gstats(int threadID,
                                ThreadRange range,
                                const std::string &binPath,
                                const std::string &metaPath)
{
    FinalBinReader reader(binPath, metaPath, g_k, range.p1Start, range.p1End);
    if(!reader.good()) return;
    bool haveCurr=false;
    uint128_t currKmer=0;
    std::vector<uint32_t> counts(g_nIDs,0);

    auto flush=[&](){
        updatePartial(g_threadPartials[threadID], counts);
    };

    while(true){
        NextRecord nr= reader.getNext();
        if(!nr.valid){
            if(haveCurr) flush();
            break;
        }
        if(!haveCurr){
            haveCurr=true;
            currKmer= nr.kmer;
            std::fill(counts.begin(), counts.end(),0);
            counts[nr.genomeID]+= nr.count;
        } else {
            if(nr.kmer==currKmer){
                counts[nr.genomeID]+= nr.count;
            } else {
                flush();
                currKmer= nr.kmer;
                std::fill(counts.begin(), counts.end(),0);
                counts[nr.genomeID]+= nr.count;
            }
        }
    }
}

static int doComputeGstats(const std::string &binPath,
                           const std::string &metaPath,
                           int numThreads)
{
    // read final.bin header
    {
        std::ifstream ifs(binPath,std::ios::binary);
        if(!ifs){
            std::cerr<<"Cannot open "<<binPath<<"\n";
            return 1;
        }
        char mg[4];
        if(!ifs.read(mg,4)){
            std::cerr<<"Cannot read magic\n";return 1;
        }
        if(std::strncmp(mg,"BINF",4)!=0){
            std::cerr<<"Not a valid BINF\n";return 1;
        }
        if(!ifs.read((char*)&g_k,sizeof(g_k))){
            std::cerr<<"Cannot read k\n";return 1;
        }
        if(g_k<1||g_k> KMER_MAX_LENGTH){
            std::cerr<<"Invalid k\n";return 1;
        }
        if(!ifs.read((char*)&g_nIDs,sizeof(g_nIDs))){
            std::cerr<<"Cannot read nIDs\n";return 1;
        }
        if(g_nIDs<1){
            std::cerr<<"No genome IDs?\n";return 1;
        }
        g_genomeNames.resize(g_nIDs);
        for(uint8_t i=0;i<g_nIDs;i++){
            uint16_t len=0;
            if(!ifs.read((char*)&len,sizeof(len))){
                std::cerr<<"Cannot read genome name length\n";return 1;
            }
            g_genomeNames[i].resize(len);
            if(!ifs.read((char*)g_genomeNames[i].data(),len)){
                std::cerr<<"Cannot read genome name\n";return 1;
            }
        }
        {
            char sentinel[4];
            if(!ifs.read(sentinel,4)){
                std::cerr<<"Cannot read ENDG\n";return 1;
            }
            if(std::strncmp(sentinel,"ENDG",4)!=0){
                std::cerr<<"Missing ENDG\n";return 1;
            }
        }
    }

    g_numThreads= std::max(numThreads,1);
    g_threadPartials.resize(g_numThreads);
    for(int t=0; t<g_numThreads; t++){
        initPartialGstats(g_threadPartials[t], g_nIDs);
    }

    int totalP1=(1<<P1_BITS);
    int chunkSize= (totalP1+g_numThreads-1)/g_numThreads;
    std::vector<ThreadRange> ranges(g_numThreads);
    for(int i=0;i<g_numThreads;i++){
        uint16_t st=(uint16_t)(i*chunkSize);
        uint16_t ed=(uint16_t)std::min<int>((i+1)*chunkSize-1, totalP1-1);
        ranges[i]={st,ed};
    }

    {
        std::vector<std::thread> threads;
        for(int i=0;i<g_numThreads;i++){
            threads.emplace_back(workerThread_gstats,i,ranges[i],binPath,metaPath);
        }
        for(auto &t: threads){
            t.join();
        }
    }

    // merge partials
    PartialGstats global;
    initPartialGstats(global, g_nIDs);
    for(int i=0;i<g_nIDs;i++){
        global.minCount[i]=(uint32_t)-1;
        global.maxCount[i]=0;
    }
    for(int t=0;t<g_numThreads;t++){
        auto &pp= g_threadPartials[t];
        for(int i=0;i<g_nIDs;i++){
            if(pp.minCount[i]< global.minCount[i]) global.minCount[i]= pp.minCount[i];
            if(pp.maxCount[i]> global.maxCount[i]) global.maxCount[i]= pp.maxCount[i];
            global.sum1[i]+= pp.sum1[i];
            global.sum2[i]+= pp.sum2[i];
            global.sum3[i]+= pp.sum3[i];
            global.nKmers[i]+=pp.nKmers[i];
        }
    }
    // merge freq
    std::vector<std::unordered_map<uint32_t,uint64_t>> globalFreq(g_nIDs);
    for(int t=0;t<g_numThreads;t++){
        auto &pp = g_threadPartials[t];
        for(int i=0;i<g_nIDs;i++){
            for(auto &kv: pp.freqMap[i]){
                globalFreq[i][kv.first]+= kv.second;
            }
        }
    }

    // final stats
    std::ofstream ofs("genome_stats.txt");
    if(!ofs){
        std::cerr<<"Cannot create genome_stats.txt\n";
        return 1;
    }
    ofs<<"#Genome\tMin\tMax\tMean\tMedian\tVariance\tSkewness\n";
    for(int i=0;i<g_nIDs;i++){
        uint64_t N= global.nKmers[i];
        uint32_t minv= (N>0? global.minCount[i] : 0);
        uint32_t maxv= global.maxCount[i];
        long double mean=0.0L;
        long double variance=0.0L;
        long double median=0.0L;
        long double skew=0.0L;
        if(N>0){
            long double s1= global.sum1[i].convert_to<long double>();
            mean= s1/(long double)N;
            long double s2= global.sum2[i].convert_to<long double>();
            long double M2= s2/(long double)N;
            variance= M2- (mean*mean);
            long double s3= global.sum3[i].convert_to<long double>();
            long double M3= s3/(long double)N;
            if(variance>0.0L){
                long double stdev= std::sqrt(variance);
                long double numer= M3 - 3.0L*mean*M2 + 2.0L*mean*mean*mean;
                long double denom= stdev*stdev*stdev;
                if(fabsl(denom)>1e-30L){
                    skew= numer/denom;
                }
            }
            // median => from freq
            uint64_t target=(N+1)/2;
            uint64_t accum=0;
            std::vector<std::pair<uint32_t,uint64_t>> fvec;
            fvec.reserve(globalFreq[i].size());
            for(auto &kv: globalFreq[i]){
                fvec.push_back(kv);
            }
            std::sort(fvec.begin(), fvec.end(),[](auto &a, auto &b){
                return a.first< b.first;
            });
            for(auto &kv: fvec){
                accum+= kv.second;
                if(accum>= target){
                    median= (long double)kv.first;
                    break;
                }
            }
        }
        std::string gName= (i< (int)g_genomeNames.size())? g_genomeNames[i]:("UnknownID#"+std::to_string(i));
        ofs<< gName<<"\t"<< minv <<"\t"<< maxv <<"\t"<< (double)mean <<"\t"<< (double)median
           <<"\t"<< (double)variance <<"\t"<< (double)skew <<"\n";
    }
    ofs.close();
    std::cout<<"Wrote genome_stats.txt.\n";
    return 0;
}

// -------------------------------------------------------------------------------------------
// --query => user provides k-mers (comma delimited or @file). We'll do a single pass to find them
// -------------------------------------------------------------------------------------------
static bool validateKmer(const std::string &kmer, int K){
    if((int)kmer.size()!=K) return false;
    for(char c: kmer){
        switch(c){
            case 'A': case 'C': case 'G': case 'T': break;
            default: return false;
        }
    }
    return true;
}

// Convert kmer string to 2-bit representation
static uint128_t encodeKmer(const std::string &kmer){
    uint128_t val=0;
    for(char c: kmer){
        val<<=2;
        switch(c){
            case 'A': val|=0; break;
            case 'C': val|=1; break;
            case 'G': val|=2; break;
            case 'T': val|=3; break;
        }
    }
    return val;
}

static int doQuery(const std::string &binPath,
                   const std::string &metaPath,
                   const std::vector<std::string> &kmersToSearch,
                   int numThreads)
{
    // read final.bin header
    {
        std::ifstream ifs(binPath,std::ios::binary);
        if(!ifs){
            std::cerr<<"Cannot open "<<binPath<<"\n";
            return 1;
        }
        char mg[4];
        if(!ifs.read(mg,4)){
            std::cerr<<"Cannot read magic\n";return 1;
        }
        if(std::strncmp(mg,"BINF",4)!=0){
            std::cerr<<"Not a valid BINF\n";return 1;
        }
        if(!ifs.read((char*)&g_k,sizeof(g_k))){
            std::cerr<<"Cannot read k\n";return 1;
        }
        if(g_k<1||g_k> KMER_MAX_LENGTH){
            std::cerr<<"Invalid k\n";return 1;
        }
        if(!ifs.read((char*)&g_nIDs,sizeof(g_nIDs))){
            std::cerr<<"Cannot read nIDs\n";return 1;
        }
        if(g_nIDs<1){
            std::cerr<<"No genome IDs\n";return 1;
        }
        g_genomeNames.resize(g_nIDs);
        for(uint8_t i=0;i<g_nIDs;i++){
            uint16_t len=0;
            if(!ifs.read((char*)&len,sizeof(len))){
                std::cerr<<"Cannot read genome name length\n";return 1;
            }
            g_genomeNames[i].resize(len);
            if(!ifs.read((char*)g_genomeNames[i].data(),len)){
                std::cerr<<"Cannot read genome name\n";return 1;
            }
        }
        {
            char sentinel[4];
            if(!ifs.read(sentinel,4)){
                std::cerr<<"Cannot read ENDG\n";return 1;
            }
            if(std::strncmp(sentinel,"ENDG",4)!=0){
                std::cerr<<"Missing ENDG\n";return 1;
            }
        }
    }

    // encode queries
    std::unordered_map<uint128_t, std::string> queryMap;
    for(auto &kmer: kmersToSearch){
        if(!validateKmer(kmer,g_k)){
            std::cerr<<"Warning: invalid k-mer '"<<kmer<<"' skipped.\n";
            continue;
        }
        uint128_t val= encodeKmer(kmer);
        queryMap[val]= kmer; 
    }
    if(queryMap.empty()){
        std::cerr<<"No valid queries.\n";
        return 0;
    }
    // We'll keep final counts in resultCounts
    std::unordered_map<uint128_t, std::vector<uint32_t>> resultCounts;
    resultCounts.reserve(queryMap.size());
    for(auto &kv: queryMap){
        resultCounts[kv.first].resize(g_nIDs,0);
    }

    // per-thread local
    struct LocalQ {
        std::unordered_map<uint128_t,std::vector<uint32_t>> localC;
    };
    std::vector<LocalQ> localVec(numThreads);
    for(int t=0;t<numThreads;t++){
        localVec[t].localC.reserve(queryMap.size());
        for(auto &kv: queryMap){
            localVec[t].localC[kv.first].resize(g_nIDs,0);
        }
    }

    auto threadFunc=[&](int threadID, ThreadRange rng){
        FinalBinReader reader(binPath, metaPath, g_k, rng.p1Start, rng.p1End);
        if(!reader.good()) return;
        bool haveCurr=false;
        uint128_t currVal=0;
        std::vector<uint32_t> counts(g_nIDs,0);

        auto flush=[&](){
            auto it= localVec[threadID].localC.find(currVal);
            if(it!= localVec[threadID].localC.end()){
                for(int i=0;i<(int)g_nIDs;i++){
                    it->second[i]+= counts[i];
                }
            }
        };

        while(true){
            NextRecord nr= reader.getNext();
            if(!nr.valid){
                if(haveCurr) flush();
                break;
            }
            if(!haveCurr){
                haveCurr=true;
                currVal= nr.kmer;
                std::fill(counts.begin(), counts.end(),0);
                counts[nr.genomeID]+= nr.count;
            } else {
                if(nr.kmer==currVal){
                    counts[nr.genomeID]+= nr.count;
                } else {
                    flush();
                    currVal= nr.kmer;
                    std::fill(counts.begin(), counts.end(),0);
                    counts[nr.genomeID]+= nr.count;
                }
            }
        }
    };

    int totalP1=(1<<P1_BITS);
    int chunkSize=(totalP1+numThreads-1)/numThreads;
    std::vector<ThreadRange> thrRanges(numThreads);
    for(int i=0;i<numThreads;i++){
        uint16_t st=(uint16_t)(i*chunkSize);
        uint16_t ed=(uint16_t)std::min<int>((i+1)*chunkSize-1, totalP1-1);
        thrRanges[i]={st,ed};
    }

    std::vector<std::thread> threads;
    for(int i=0;i<numThreads;i++){
        threads.emplace_back(threadFunc, i, thrRanges[i]);
    }
    for(auto &t: threads){
        t.join();
    }

    // merge local => global
    for(int t=0;t<numThreads;t++){
        for(auto &kv: localVec[t].localC){
            auto &glob= resultCounts[kv.first];
            auto &loc = kv.second;
            for(int i=0;i<(int)g_nIDs;i++){
                glob[i]+= loc[i];
            }
        }
    }

    // print
    for(auto &kv: queryMap){
        auto &kVal = kv.first;
        auto &kStr = kv.second;
        auto &cnts = resultCounts[kVal];
        std::cout<< kStr <<": ";
        bool first=true;
        for(int i=0;i<(int)g_nIDs;i++){
            if(!first) std::cout<<", ";
            first=false;
            if(i<(int)g_genomeNames.size()){
                std::cout<< g_genomeNames[i] <<": "<< cnts[i];
            } else {
                std::cout<<"UnknownID#"<<i <<": "<< cnts[i];
            }
        }
        std::cout<<"\n";
    }
    return 0;
}

// -------------------------------------------------------------------------------------------
// --query_regex => decode each k-mer, check if it matches the user-provided regex, collect
// -------------------------------------------------------------------------------------------
struct MatchedKmer {
    std::string kmer;
    std::vector<uint32_t> counts;
};

static void workerThread_regex(int threadID,
                               ThreadRange range,
                               const std::string &binPath,
                               const std::string &metaPath,
                               const std::regex &re,
                               std::vector<MatchedKmer> &outVec)
{
    FinalBinReader reader(binPath, metaPath, g_k, range.p1Start, range.p1End);
    if(!reader.good()) return;

    bool haveCurr=false;
    uint128_t currVal=0;
    std::vector<uint32_t> counts(g_nIDs,0);

    auto flush=[&](){
        // decode
        KmerResult kr = decodeKmer128(currVal,g_k);
        // check regex
        if(std::regex_match(kr.kmer, re)){
            // store
            MatchedKmer mk;
            mk.kmer= kr.kmer;
            mk.counts= counts;
            outVec.push_back(std::move(mk));
        }
    };

    while(true){
        NextRecord nr= reader.getNext();
        if(!nr.valid){
            if(haveCurr) flush();
            break;
        }
        if(!haveCurr){
            haveCurr=true;
            currVal= nr.kmer;
            std::fill(counts.begin(), counts.end(),0);
            counts[nr.genomeID]+= nr.count;
        } else {
            if(nr.kmer==currVal){
                counts[nr.genomeID]+= nr.count;
            } else {
                flush();
                currVal= nr.kmer;
                std::fill(counts.begin(),counts.end(),0);
                counts[nr.genomeID]+= nr.count;
            }
        }
    }
}

static int doQueryRegex(const std::string &binPath,
                        const std::string &metaPath,
                        int threads,
                        const std::string &pattern)
{
    // read final.bin header
    {
        std::ifstream ifs(binPath,std::ios::binary);
        if(!ifs){
            std::cerr<<"Cannot open "<<binPath<<"\n";
            return 1;
        }
        char mg[4];
        if(!ifs.read(mg,4)){
            std::cerr<<"Cannot read magic\n";return 1;
        }
        if(std::strncmp(mg,"BINF",4)!=0){
            std::cerr<<"Not a valid BINF\n";return 1;
        }
        if(!ifs.read((char*)&g_k,sizeof(g_k))){
            std::cerr<<"Cannot read k\n";return 1;
        }
        if(g_k<1||g_k> KMER_MAX_LENGTH){
            std::cerr<<"Invalid k\n";return 1;
        }
        if(!ifs.read((char*)&g_nIDs,sizeof(g_nIDs))){
            std::cerr<<"Cannot read nIDs\n";return 1;
        }
        if(g_nIDs<1){
            std::cerr<<"No genome IDs\n";return 1;
        }
        g_genomeNames.resize(g_nIDs);
        for(uint8_t i=0;i<g_nIDs;i++){
            uint16_t len=0;
            if(!ifs.read((char*)&len,sizeof(len))){
                std::cerr<<"Cannot read genome name length\n";return 1;
            }
            g_genomeNames[i].resize(len);
            if(!ifs.read((char*)g_genomeNames[i].data(),len)){
                std::cerr<<"Cannot read genome name\n";return 1;
            }
        }
        {
            char sentinel[4];
            if(!ifs.read(sentinel,4)){
                std::cerr<<"Cannot read ENDG\n";return 1;
            }
            if(std::strncmp(sentinel,"ENDG",4)!=0){
                std::cerr<<"Missing ENDG\n";return 1;
            }
        }
    }

    // compile regex
    std::regex re(pattern);

    // We'll do a multi-thread approach
    int totalP1=(1<<P1_BITS);
    int chunkSize=(totalP1+threads-1)/threads;
    std::vector<ThreadRange> ranges(threads);
    for(int i=0;i<threads;i++){
        uint16_t st=(uint16_t)(i*chunkSize);
        uint16_t ed=(uint16_t)std::min<int>((i+1)*chunkSize-1, totalP1-1);
        ranges[i]={st,ed};
    }

    // We'll store each thread's matched results in a separate vector, then combine
    std::vector<std::vector<MatchedKmer>> allThreadResults(threads);

    std::vector<std::thread> tvec;
    for(int i=0;i<threads;i++){
        tvec.emplace_back(workerThread_regex,
                          i,
                          ranges[i],
                          binPath,
                          metaPath,
                          std::cref(re),
                          std::ref(allThreadResults[i]));
    }
    for(auto &th: tvec){
        th.join();
    }

    // combine
    // For simplicity, we won't reorder them. We'll just print them in the order found.
    for(int i=0;i<threads;i++){
        for(auto &mk: allThreadResults[i]){
            // format: "kmer: G1: c1, G2: c2"
            std::cout << mk.kmer << ": ";
            bool first=true;
            for(int j=0;j<(int)g_nIDs;j++){
                uint32_t c= mk.counts[j];
                if(!first) std::cout<<", ";
                first=false;
                std::string nm = (j<(int)g_genomeNames.size()? g_genomeNames[j]:("UnknownID#"+std::to_string(j)));
                std::cout<< nm <<": "<< c;
            }
            std::cout<<"\n";
        }
    }
    return 0;
}

// -------------------------------------------------------------------------------------------
// MAIN
//   usage is now:
//   maf_counter_tools --std <topCount> --threads <N> --metadata_file <metaFile> --binary_database <binFile>
//   maf_counter_tools --expr <expression> --threads <N> --metadata_file <metaFile> --binary_database <binFile>
//   maf_counter_tools --gstats --threads <N> --metadata_file <metaFile> --binary_database <binFile>
//   maf_counter_tools --query <kmersOr@file> --threads <N> --metadata_file <metaFile> --binary_database <binFile>
//   maf_counter_tools --query_regex <REGEX> --threads <N> --metadata_file <metaFile> --binary_database <binFile>
// -------------------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
    if(argc<2){
        std::cerr<<"Usage:\n"
                 <<"  "<<argv[0]<<" --std <topCount> --threads <N> --metadata_file <META> --binary_database <BIN>\n"
                 <<"  "<<argv[0]<<" --expr <EXPR> --threads <N> --metadata_file <META> --binary_database <BIN>\n"
                 <<"  "<<argv[0]<<" --gstats --threads <N> --metadata_file <META> --binary_database <BIN>\n"
                 <<"  "<<argv[0]<<" --query <kmerListOr@file> --threads <N> --metadata_file <META> --binary_database <BIN>\n"
                 <<"  "<<argv[0]<<" --query_regex <REGEX> --threads <N> --metadata_file <META> --binary_database <BIN>\n"
                 <<"\n"
                 <<"Examples:\n"
                 <<"  "<<argv[0]<<" --std 20 --threads 4 --metadata_file final.metadata --binary_database final.bin\n"
                 <<"  "<<argv[0]<<" --expr \"CHM13>100 && GRCh38>20\" --threads 6 --metadata_file final.metadata --binary_database final.bin\n"
                 <<"  "<<argv[0]<<" --gstats --threads 4 --metadata_file final.metadata --binary_database final.bin\n"
                 <<"  "<<argv[0]<<" --query \"ACGTACGT,AAAAAAAC\" --threads 2 --metadata_file final.metadata --binary_database final.bin\n"
                 <<"  "<<argv[0]<<" --query_regex \"([gG]{3,}\\w{1,7}){3,}[gG]{3,}\" --threads 2 --metadata_file final.metadata --binary_database final.bin\n";
        return 1;
    }

    // We'll parse arguments more manually
    std::string mode;
    int topCount=-1;
    std::string exprStr;
    std::string queryStr;
    std::string regexStr;
    int numThreads=-1;
    std::string metaFile = "final.metadata";
    std::string binFile = "final.bin";

    

    // We have possible flags: 
    //  --std <val>, --expr <val>, --gstats, --query <val>, --query_regex <val>
    //  --threads <val>, --metadata_file <val>, --binary_database <val>
    // We'll parse them in order
    for(int i=1;i<argc;i++){
        std::string arg= argv[i];
        if(arg=="--std"){
            mode="std";
            if(i+1<argc){
                topCount= std::stoi(argv[++i]);
            } else {
                std::cerr<<"Missing argument for --std\n";
                return 1;
            }
        }
        else if(arg=="--expr"){
            mode="expr";
            if(i+1<argc){
                exprStr= argv[++i];
            } else {
                std::cerr<<"Missing argument for --expr\n";
                return 1;
            }
        }
        else if(arg=="--gstats"){
            mode="gstats";
        }
        else if(arg=="--query"){
            mode="query";
            if(i+1<argc){
                queryStr= argv[++i];
            } else {
                std::cerr<<"Missing argument for --query\n";
                return 1;
            }
        }
        else if(arg=="--query_regex"){
            mode="query_regex";
            if(i+1<argc){
                regexStr= argv[++i];
            } else {
                std::cerr<<"Missing argument for --query_regex\n";
                return 1;
            }
        }
        else if(arg=="--threads"){
            if(i+1<argc){
                numThreads= std::stoi(argv[++i]);
            } else {
                std::cerr<<"Missing argument for --threads\n";
                return 1;
            }
        }
        else if(arg=="--metadata_file"){
            if(i+1<argc){
                metaFile= argv[++i];
            } else {
                std::cerr<<"Missing argument for --metadata_file\n";
                return 1;
            }
        }
        else if(arg=="--binary_database"){
            if(i+1<argc){
                binFile= argv[++i];
            } else {
                std::cerr<<"Missing argument for --binary_database\n";
                return 1;
            }
        }
        else {
            std::cerr<<"Unknown argument: "<<arg<<"\n";
            return 1;
        }
    }

    if(mode.empty()){
        std::cerr<<"You must specify one of --std / --expr / --gstats / --query / --query_regex\n";
        return 1;
    }
    if(numThreads<1){
        std::cerr<<"You must specify --threads <val> with val>0\n";
        return 1;
    }
    if(metaFile.empty()){
        std::cerr<<"You must specify --metadata_file <META>\n";
        return 1;
    }
    if(binFile.empty()){
        std::cerr<<"You must specify --binary_database <BIN>\n";
        return 1;
    }

    // Now dispatch
    if(mode=="std"){
        if(topCount<1){
            std::cerr<<"Invalid topCount for --std\n";
            return 1;
        }
        return doComputeTopStds(binFile, metaFile, topCount, numThreads);
    }
    else if(mode=="expr"){
        if(exprStr.empty()){
            std::cerr<<"No expression for --expr\n";
            return 1;
        }
        return doComputeExpr(binFile, metaFile, numThreads, exprStr);
    }
    else if(mode=="gstats"){
        return doComputeGstats(binFile, metaFile, numThreads);
    }
    else if(mode=="query"){
        if(queryStr.empty()){
            std::cerr<<"No argument for --query\n";
            return 1;
        }
        // parse queryStr: if starts with '@', read from file, else comma-split
        std::vector<std::string> kms;
        if(!queryStr.empty() && queryStr[0]=='@'){
            std::string f= queryStr.substr(1);
            std::ifstream qifs(f);
            if(!qifs){
                std::cerr<<"Cannot open k-mer query file "<<f<<"\n";
                return 1;
            }
            std::string line;
            while(std::getline(qifs,line)){
                if(!line.empty()) kms.push_back(line);
            }
            qifs.close();
        } else {
            std::stringstream ss(queryStr);
            std::string tok;
            while(std::getline(ss,tok,',')){
                if(!tok.empty()) kms.push_back(tok);
            }
        }
        return doQuery(binFile, metaFile, kms, numThreads);
    }
    else if(mode=="query_regex"){
        if(regexStr.empty()){
            std::cerr<<"No argument for --query_regex\n";
            return 1;
        }
        return doQueryRegex(binFile, metaFile, numThreads, regexStr);
    }
    // fallback
    std::cerr<<"Unrecognized mode.\n";
    return 1;
}
