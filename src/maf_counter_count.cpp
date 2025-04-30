#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <atomic>
#include <unordered_map>
#include <sstream>
#include <algorithm>
#include <random>
#include <cstdint>
#include <cstdlib>
#include <chrono>
#include <cstring>
#include <filesystem>
#include <regex>

// google dense hash
#include <google/dense_hash_map>
#include <unordered_set>
#include <limits>
#include <zlib.h>
#include <sys/mman.h>
#include <fcntl.h>      // for open, O_*
#include <sys/sendfile.h> // for sendfile()
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

// ----  NEW: Boost Multiprecision for 128-bit support  ----
#include <boost/multiprecision/cpp_int.hpp>
// --------------------------------------------------------
namespace bmp = boost::multiprecision;
using uint128_t = bmp::uint128_t;

// Global mutex for printing timing results.
std::mutex coutMutex;

// ---------------------------------------------
// Common definitions
// ---------------------------------------------
static const int   P1_BITS            = 10;
static const int   P2_BITS            = 10;
static const size_t BIN_CAPACITY      = 1ULL << 14;
static const size_t PACKAGE_THRESHOLD = 500000000;

static const std::unordered_set<uint16_t> outlierBins = {
    1023,0,3,255,1020,831,12,768,128,1015,15,1011,192,63,1008,1022,
    1021,512,960,975,48,32,8,991,256,895,2,511,14,319,1019,959,
    4,64,1012,896,252,771,207,51,136,887,338,819,204,1,767,490,
    983,60,160,783
};
static std::mutex outlierMutex;

// Expanded maximum k-mer length to 64
static const int   KMER_MAX_LENGTH  = 64;  

// New global variables for q-line filtering
static int g_minQVal = -1;  
static int g_maxQVal = -1;  
std::vector<uint64_t> g_binFileSizes(1 << P1_BITS, 0);
static double g_minAScore = -std::numeric_limits<double>::infinity();
static double g_maxAScore =  std::numeric_limits<double>::infinity();
static std::unordered_set<std::string> g_includedGenomeNames;
static bool g_useGenomeIdFilter = false;
static std::mutex g_genomeMapMutex;
static std::unordered_map<std::string, uint8_t> g_genomeMap;
static std::vector<std::string> g_genomeNames;  
static uint8_t g_nextGenomeID = 0;
static int g_k = 0;
static int g_numReaders = 1;
static int g_numPMs = 1;
static bool g_purgeIntermediate = true; 
static std::string g_binaryOutputFile = "final.bin";
static bool g_useCanonical = false;  

// New global variables for directories
static std::string g_tempFilesDir = ".";   // Where intermediate bin files are stored
static std::string g_outputDirectory = "."; // Where final.bin and final.metadata are written

struct TimeStats {
    double readTime{};
    double sortTime{};
    double writeTime{};
};

// Convert a single-character quality into an integer.
static int parseQVal(char c)
{
    if (c >= '0' && c <= '9')
        return c - '0';
    if (c == 'F' || c == 'f')
        return 15;
    return -1;
}

// Check if a single-character quality is in the allowed range.
static bool isQValid(char c)
{
    if (g_minQVal < 0 && g_maxQVal < 0) {
        return true;
    }
    int val = parseQVal(c);
    if (val < 0) {
        return false;
    }
    if (g_minQVal >= 0 && val < g_minQVal) {
        return false;
    }
    if (g_maxQVal >= 0 && val > g_maxQVal) {
        return false;
    }
    return true;
}

// ---------------------------------------------
// .gz handling logic
// ---------------------------------------------
static std::string decompressGZ(const std::string &inFile)
{
    std::string tmpFileName = "temp_decompressed_";
    tmpFileName += std::to_string(std::chrono::steady_clock::now().time_since_epoch().count());
    tmpFileName += ".maf";

    gzFile in = gzopen(inFile.c_str(), "rb");
    if (!in) {
        std::cerr << "Error: could not open compressed file " << inFile << "\n";
        return {};
    }

    FILE *out = std::fopen(tmpFileName.c_str(), "wb");
    if (!out) {
        std::cerr << "Error: could not create temporary file " << tmpFileName << "\n";
        gzclose(in);
        return {};
    }

    const size_t BUF_SIZE = 65536;
    std::vector<char> buffer(BUF_SIZE);
    int bytesRead = 0;
    while ((bytesRead = gzread(in, buffer.data(), static_cast<unsigned int>(buffer.size()))) > 0) {
        std::fwrite(buffer.data(), 1, bytesRead, out);
    }

    gzclose(in);
    std::fclose(out);

    return tmpFileName;
}

static bool endsWithGZ(const std::string &filename)
{
    if (filename.size() < 3) return false;
    return (filename.compare(filename.size() - 3, 3, ".gz") == 0);
}

static void removeFileIfExists(const std::string &filename)
{
    std::error_code ec;
    std::filesystem::remove(filename, ec);
}

// ---------------------------------------------
// Common k-mer ID retrieval and other routines
// ---------------------------------------------
static uint8_t getGenomeID(const std::string &genomeName)
{
    std::lock_guard<std::mutex> lock(g_genomeMapMutex);
    auto it = g_genomeMap.find(genomeName);
    if (it != g_genomeMap.end()) {
        return it->second;
    } else {
        if (g_nextGenomeID == 255) {
            std::cerr << "Too many distinct genomes\n";
            std::exit(1);
        }
        uint8_t newID = g_nextGenomeID++;
        g_genomeMap[genomeName] = newID;
        if (newID >= g_genomeNames.size()) {
            g_genomeNames.resize(newID + 1);
        }
        g_genomeNames[newID] = genomeName;
        return newID;
    }
}

static void findMAFChunkBoundaries(const std::string &filename,
                                   int numChunks,
                                   std::vector<std::streampos> &chunkStartPositions)
{
    std::ifstream mafFile(filename);
    if (!mafFile) {
        std::cerr << "Error opening MAF file\n";
        std::exit(1);
    }
    mafFile.seekg(0, std::ios::end);
    std::streampos fileSize = mafFile.tellg();
    mafFile.seekg(0, std::ios::beg);
    std::vector<std::streampos> approx(numChunks + 1);
    for (int i = 0; i <= numChunks; i++) {
        approx[i] = fileSize * i / numChunks;
    }
    chunkStartPositions.resize(numChunks + 1);
    chunkStartPositions[0] = 0;
    for (int i = 1; i < numChunks; i++) {
        mafFile.seekg(approx[i]);
        std::string line;
        while (mafFile.tellg() < fileSize && std::getline(mafFile, line)) {
            if (line.empty()) continue;
            if (line[0] == 'a') {
                std::streampos currentPos = mafFile.tellg();
                std::streamoff lineLen = (std::streamoff)line.size() + 1;
                chunkStartPositions[i] = currentPos - lineLen;
                break;
            }
        }
    }
    chunkStartPositions[numChunks] = fileSize;
}

static std::string removeDashes(const std::string &seq)
{
    std::string result;
    result.reserve(seq.size());
    for (char c : seq) {
        if (c != '-') result.push_back(c);
    }
    return result;
}

// ---- decodeKmer for 128 bits ----
#include <boost/multiprecision/cpp_int.hpp> // repeated but for clarity
static std::string decodeKmer(uint128_t kmerVal, int K)
{
    std::string result(K, 'N');
    for (int i = 0; i < K; i++) {
        int shift = 2 * (K - 1 - i);
        uint128_t code = (kmerVal >> shift) & 3;
        switch (code.convert_to<unsigned long>()) {
            case 0: result[i] = 'A'; break;
            case 1: result[i] = 'C'; break;
            case 2: result[i] = 'G'; break;
            case 3: result[i] = 'T'; break;
        }
    }
    return result;
}

// ---- Reverse complement logic (128-bit) ----
static uint128_t reverseComplement128(uint128_t kmer, int k)
{
    uint128_t rc = 0;
    for (int i = 0; i < k; i++) {
        uint128_t bits = (kmer & 3U); 
        bits ^= 3U;
        rc = (rc << 2) | bits;
        kmer >>= 2;
    }
    return rc;
}

static uint128_t canonicalKmer128(uint128_t kmer, int k)
{
    uint128_t rc = reverseComplement128(kmer, k);
    return (rc < kmer) ? rc : kmer;
}

// ---------------------------------------------
// K <= 10 approach (32-bit k-mers)
// ---------------------------------------------
struct LocalHashData {
    google::dense_hash_map<uint32_t, std::array<uint32_t, 256>> map;
};

struct SmallKReaderArgs {
    std::string filename;
    std::streampos startPos;
    std::streampos endPos;
    int k;
    LocalHashData* localData;
    std::atomic<uint64_t>* totalKmers;
    int readerId;
};

static int parseQVal(char c);  // (already defined above)

static void smallKReaderFunc(SmallKReaderArgs args)
{
    auto start = std::chrono::steady_clock::now();

    args.localData->map.set_empty_key(0xFFFFFFFF);
    args.localData->map.set_deleted_key(0xFFFFFFFE);

    std::ifstream mafFile(args.filename);
    mafFile.seekg(args.startPos);

    uint64_t localCount = 0;
    std::string line;
    double currentBlockScore = 0.0;

    while (true) {
        if (mafFile.tellg() >= args.endPos) break;
        if (!std::getline(mafFile, line)) break;
        if (line.empty()) continue;
        
        if (line[0] == 'a'){
            std::size_t pos = line.find("score=");
            if (pos != std::string::npos) {
                std::string scoreStr = line.substr(pos + 6);
                try {
                    currentBlockScore = std::stod(scoreStr);
                } catch (...) {
                    currentBlockScore = 0.0; 
                }
            } else {
                currentBlockScore = 0.0;
            }
            continue;
        }

        if (line[0] != 's') {
            continue;
        }

        if (currentBlockScore < g_minAScore || currentBlockScore > g_maxAScore) {
            continue;
        }

        std::istringstream iss(line);
        std::string token, src, seq;
        iss >> token >> src; 
        for (int i = 0; i < 4; i++) iss >> token;
        iss >> seq;

        size_t dotPos = src.find('.');
        std::string genomeOnly = (dotPos == std::string::npos) ? src : src.substr(0, dotPos);
        if (g_useGenomeIdFilter && g_includedGenomeNames.find(genomeOnly) == g_includedGenomeNames.end()) {
            continue;
        }
        uint8_t gID = getGenomeID(genomeOnly);

        std::string cleanedQ;
        if ((g_minQVal >= 0 || g_maxQVal >= 0)) {
            std::streampos savePos = mafFile.tellg();
            std::string possibleQ;
            if (std::getline(mafFile, possibleQ)) {
                if (!possibleQ.empty() && possibleQ[0] == 'q') {
                    std::istringstream qiss(possibleQ);
                    std::string qtoken, qsrc, qseq;
                    qiss >> qtoken >> qsrc;
                    for (int i = 0; i < 4; i++) qiss >> qtoken;
                    qiss >> qseq;
                    cleanedQ = removeDashes(qseq);
                    std::transform(cleanedQ.begin(), cleanedQ.end(), cleanedQ.begin(),
               [](char c){ return std::toupper(c); });
                } else {
                    mafFile.seekg(savePos);
                }
            }
        }

        std::string cleaned = removeDashes(seq);
        std::transform(cleaned.begin(), cleaned.end(), cleaned.begin(),
               [](char c){ return std::toupper(c); });
        
        if (cleaned.size() < (size_t)args.k) continue;

        for (size_t pos = 0; pos + args.k <= cleaned.size(); pos++) {
            uint32_t kmer32 = 0;
            bool valid = true;
            for (int i = 0; i < args.k; i++) {
                char c = cleaned[pos + i];
                uint32_t val;
                switch (c) {
                    case 'A': case 'a': val = 0; break;
                    case 'C': case 'c': val = 1; break;
                    case 'G': case 'g': val = 2; break;
                    case 'T': case 't': val = 3; break;
                    default:  valid = false; break;
                }
                if (!valid) break;
                if (!cleanedQ.empty()) {
                    if (!isQValid(cleanedQ[pos + i])) {
                        valid = false;
                        break;
                    }
                }
                kmer32 = (kmer32 << 2) | val;
            }
            if (!valid) continue;

            if (g_useCanonical) {
                uint32_t rc = 0;
                uint32_t tmp = kmer32;
                for (int i = 0; i < args.k; i++) {
                    uint32_t bits = tmp & 3U;
                    bits ^= 3U;
                    rc = (rc << 2) | bits;
                    tmp >>= 2;
                }
                if (rc < kmer32) {
                    kmer32 = rc;
                }
            }

            auto it = args.localData->map.find(kmer32);
            if (it == args.localData->map.end()) {
                std::array<uint32_t, 256> arr{};
                arr[gID]++;
                args.localData->map.insert({kmer32, arr});
            } else {
                it->second[gID]++;
            }
            localCount++;
        }
    }
    args.totalKmers->fetch_add(localCount, std::memory_order_relaxed);

    auto end = std::chrono::steady_clock::now();
    double duration = std::chrono::duration<double>(end - start).count();
    {
        std::lock_guard<std::mutex> lock(coutMutex);
        std::cout << "SmallKReader thread " << args.readerId << " finished in " << duration << " seconds.\n";
    }
}

// ---------------------------------------------
// K > 10 aggregator approach (128-bit kmers)
// ---------------------------------------------
struct KmerRecord {
    uint128_t kmer;
    uint8_t   genomeID;
};

static inline uint16_t extractP1(const uint128_t &kmer, int k)
{
    int totalBits = 2*k;
    int shift = totalBits - P1_BITS;
    uint128_t top = (kmer >> shift) & 0x3FFU; 
    return (uint16_t)top.convert_to<uint16_t>();
}

static inline uint16_t extractP2(const uint128_t &kmer, int k)
{
    int totalBits = 2*k;
    int shiftAfterP1 = totalBits - P1_BITS;
    int shiftP2      = shiftAfterP1 - P2_BITS;
    uint128_t val = (kmer >> shiftP2) & 0x3FFU;
    return (uint16_t)val.convert_to<uint16_t>();
}

struct BinBuffer {
    uint16_t p1_id;
    std::vector<KmerRecord> records;
};

struct SuffixGroup {
    uint16_t p2_val;
    uint8_t suffixMode;
    std::vector<uint32_t> mode0;
    std::vector<uint64_t> mode1;
    std::vector<uint128_t> mode3;
};

struct BinChunk {
    uint16_t p1_id;
    std::vector<SuffixGroup> groups;
};

// LSD-based sorts for 32,64,128
static inline void LSD32(std::vector<uint32_t>& a) {
    if (a.empty()) return;
    constexpr int P=4, R=256;
    std::vector<uint32_t> tmp(a.size());
    for (int p=0; p<P; ++p) {
        uint32_t cnt[R] = {0}, pos[R];
        for (auto v : a) ++cnt[(v >> (p*8)) & 0xFF];
        pos[0] = 0;
        for (int i=1; i<R; ++i) pos[i] = pos[i-1] + cnt[i-1];
        for (auto v : a) tmp[pos[(v >> (p*8)) & 0xFF]++] = v;
        a.swap(tmp);
    }
}
static inline void LSD64(std::vector<uint64_t>& a) {
    if (a.empty()) return;
    constexpr int P=8, R=256;
    std::vector<uint64_t> tmp(a.size());
    for (int p=0; p<P; ++p) {
        uint64_t cnt[R] = {0}, pos[R];
        for (auto v : a) ++cnt[(v >> (p*8)) & 0xFF];
        pos[0] = 0;
        for (int i=1; i<R; ++i) pos[i] = pos[i-1] + cnt[i-1];
        for (auto v : a) tmp[pos[(v >> (p*8)) & 0xFF]++] = v;
        a.swap(tmp);
    }
}
static inline void LSD128(std::vector<uint128_t>& a) {
    if (a.empty()) return;
    constexpr int P=16, R=256;
    std::vector<uint128_t> tmp(a.size());
    for (int p=0; p<P; ++p) {
        uint64_t cnt[R] = {0}, pos[R];
        for (auto &v : a) ++cnt[static_cast<uint8_t>((v >> (p*8)) & 0xFF)];
        pos[0] = 0;
        for (int i=1; i<R; ++i) pos[i] = pos[i-1] + cnt[i-1];
        for (auto &v : a)
            tmp[pos[static_cast<uint8_t>((v >> (p*8)) & 0xFF)]++] = v;
        a.swap(tmp);
    }
}

static uint8_t determineSuffixMode(int remainder_bits)
{
    if (remainder_bits <= 32) return 0;
    if (remainder_bits <= 64) return 1;
    return 3;
}

static BinChunk partialSortAndCompact(const BinBuffer &binBuf, int k)
{
    BinChunk chunk;
    chunk.p1_id = binBuf.p1_id;
    static const size_t P2_BUCKETS = 1ULL << P2_BITS;

    std::vector<std::vector<KmerRecord>> buckets(P2_BUCKETS);

    int remainder_bits = 2*k - P1_BITS - P2_BITS; 
    uint8_t suffixMode = determineSuffixMode(remainder_bits);

    for (auto &rec : binBuf.records) {
        uint16_t p2 = extractP2(rec.kmer, k);
        buckets[p2].push_back(rec);
    }

    chunk.groups.reserve(P2_BUCKETS);
    for (size_t i = 0; i < P2_BUCKETS; i++) {
        if (buckets[i].empty()) continue;
        SuffixGroup grp;
        grp.p2_val    = (uint16_t)i;
        grp.suffixMode = suffixMode;

        switch (suffixMode) {
            case 0: grp.mode0.reserve(buckets[i].size()); break;
            case 1: grp.mode1.reserve(buckets[i].size()); break;
            case 3: grp.mode3.reserve(buckets[i].size()); break;
        }

        for (auto &rec : buckets[i]) {
            uint128_t mask = ((uint128_t)1 << remainder_bits) - 1;
            uint128_t rem128 = rec.kmer & mask;
            switch (suffixMode) {
                case 0: {
                    uint32_t rem32 = rem128.convert_to<uint32_t>();
                    uint32_t val32 = (rem32 << 8) | rec.genomeID;
                    grp.mode0.push_back(val32);
                } break;
                case 1: {
                    uint64_t rem64 = rem128.convert_to<uint64_t>();
                    uint64_t val64 = (rem64 << 8) | rec.genomeID;
                    grp.mode1.push_back(val64);
                } break;
                case 3: {
                    uint128_t val128 = (rem128 << 8) | rec.genomeID;
                    grp.mode3.push_back(val128);
                } break;
            }
        }
        chunk.groups.push_back(std::move(grp));
    }

    return chunk;
}

// -------------------------------------
// Package Manager
// -------------------------------------
struct PMBinAccum {
    std::vector<SuffixGroup> groups;
    size_t totalCount = 0;
};

struct PackageManager {
    std::mutex pmMutex;
    std::condition_variable pmCV;
    std::queue<BinChunk> pmQueue;
    bool pmReadersDone = false;
    std::unordered_map<uint16_t, PMBinAccum> pmData;
};

static void mergeBinChunkIntoPackage(BinChunk &chunk,
                                     std::unordered_map<uint16_t, PMBinAccum> &pmData)
{
    auto &acc = pmData[chunk.p1_id];
    for (auto &grp : chunk.groups) {
        auto it = std::lower_bound(
            acc.groups.begin(), acc.groups.end(), grp.p2_val,
            [](const SuffixGroup &s, uint16_t val){return s.p2_val < val;}
        );
        if (it != acc.groups.end() && it->p2_val == grp.p2_val) {
            switch (grp.suffixMode) {
                case 0:
                    it->mode0.insert(it->mode0.end(),
                                     std::make_move_iterator(grp.mode0.begin()),
                                     std::make_move_iterator(grp.mode0.end()));
                    acc.totalCount += grp.mode0.size();
                    break;
                case 1:
                    it->mode1.insert(it->mode1.end(),
                                     std::make_move_iterator(grp.mode1.begin()),
                                     std::make_move_iterator(grp.mode1.end()));
                    acc.totalCount += grp.mode1.size();
                    break;
                case 3:
                    it->mode3.insert(it->mode3.end(),
                                     std::make_move_iterator(grp.mode3.begin()),
                                     std::make_move_iterator(grp.mode3.end()));
                    acc.totalCount += grp.mode3.size();
                    break;
            }
        } else {
            it = acc.groups.insert(it, std::move(grp));
            switch (it->suffixMode) {
                case 0: acc.totalCount += it->mode0.size(); break;
                case 1: acc.totalCount += it->mode1.size(); break;
                case 3: acc.totalCount += it->mode3.size(); break;
            }
        }
    }
}

static void writePackageToDisk(uint16_t p1_id, PMBinAccum &accum, int k)
{
    std::uniform_int_distribution<int> dist(100000,999999);
    static thread_local std::mt19937 rng(std::random_device{}());
    int rnd = dist(rng);
    std::string fname = g_tempFilesDir + "/bin_" + std::to_string(p1_id) + "_" + std::to_string(rnd) + ".bin";

    std::ofstream ofs(fname, std::ios::binary);
    if(!ofs) {
        std::cerr << "Error creating " << fname << "\n";
        return;
    }
    struct FileHeader {
        char     magic[4];
        uint32_t k;
        uint16_t p1_bits;
        uint16_t p2_bits;
        uint8_t  genome_id_bits;
        uint8_t  suffix_mode; 
        uint8_t  reserved[2];
    } hdr;

    std::memcpy(hdr.magic, "KMHD", 4);
    hdr.k            = (uint32_t)k;
    hdr.p1_bits      = P1_BITS;
    hdr.p2_bits      = P2_BITS;
    hdr.genome_id_bits = 8;

    uint8_t anyMode = 0;
    if (!accum.groups.empty()) {
        anyMode = accum.groups.front().suffixMode;
    }
    hdr.suffix_mode  = anyMode;
    hdr.reserved[0]  = 0; 
    hdr.reserved[1]  = 0;
    ofs.write((const char*)&hdr, sizeof(hdr));

    const char * sentinel1 = "kmhd";
    ofs.write(sentinel1, 4);
    uint32_t groupCount = (uint32_t)accum.groups.size();
    ofs.write((const char*)&groupCount, sizeof(groupCount));
    for (auto &g : accum.groups) {
        ofs.write((const char*)&g.p2_val, sizeof(g.p2_val));
        switch (g.suffixMode) {
            case 0: {
                uint32_t sz = (uint32_t)g.mode0.size();
                ofs.write((const char*)&sz, sizeof(sz));
            } break;
            case 1: {
                uint32_t sz = (uint32_t)g.mode1.size();
                ofs.write((const char*)&sz, sizeof(sz));
            } break;
            case 3: {
                uint32_t sz = (uint32_t)g.mode3.size();
                ofs.write((const char*)&sz, sizeof(sz));
            } break;
        }
    }
    const char* sentinel2 = "kmhd";
    ofs.write(sentinel2, 4);
    for (auto &g : accum.groups) {
        switch (g.suffixMode) {
            case 0:
                ofs.write((const char*)g.mode0.data(), g.mode0.size() * sizeof(uint32_t));
                break;
            case 1:
                for (auto &val : g.mode1) {
                    ofs.write((const char*)&val, sizeof(uint64_t));
                }
                break;
            case 3:
                for (auto &val : g.mode3) {
                    unsigned char bytes[16];
                    for (int i = 0; i < 16; i++) {
                        uint64_t part = (val >> (8ULL * i)).convert_to<uint64_t>();
                        bytes[i] = static_cast<unsigned char>(part & 0xFF);
                    }
                    ofs.write((char*)bytes, 16);
                }
                break;
        }
    }
    ofs.close();
    accum.groups.clear();
    accum.totalCount = 0;
}

static void PMThreadLoop(PackageManager &pm, int k, int pmId)
{
    auto start = std::chrono::steady_clock::now();
    double pmWriteTime = 0.0;  // total time spent in writePackageToDisk()
    size_t BIN_THRESHOLD = (PACKAGE_THRESHOLD / (1 << P1_BITS)) * 5;
    while (true) {
        std::unique_lock<std::mutex> lock(pm.pmMutex);
        pm.pmCV.wait(lock, [&pm]{return !pm.pmQueue.empty() || pm.pmReadersDone;});
        while (!pm.pmQueue.empty()) {
            BinChunk c = std::move(pm.pmQueue.front());
            pm.pmQueue.pop();
            lock.unlock();

            mergeBinChunkIntoPackage(c, pm.pmData);
            auto &acc = pm.pmData[c.p1_id];
            if (acc.totalCount >= BIN_THRESHOLD) {
                 auto w0 = std::chrono::steady_clock::now();
                writePackageToDisk(c.p1_id, acc, k);
               auto w1 = std::chrono::steady_clock::now();
                pmWriteTime += std::chrono::duration<double>(w1 - w0).count();
            }
            lock.lock();
        }
        if (pm.pmReadersDone && pm.pmQueue.empty()) break;
    }
    for (auto &kv : pm.pmData) {
        if (kv.second.totalCount > 0) {
             auto w0 = std::chrono::steady_clock::now();
            writePackageToDisk(kv.first, kv.second, k);
            auto w1 = std::chrono::steady_clock::now();
           pmWriteTime += std::chrono::duration<double>(w1 - w0).count();
        }
    }
        auto end = std::chrono::steady_clock::now();
        double totalDuration = std::chrono::duration<double>(end - start).count();
        {
            std::lock_guard<std::mutex> lock(coutMutex);
            std::cout << "PM thread " << pmId
                      << " total run time: " << totalDuration << " s, "
                      << "time in writePackageToDisk(): " << pmWriteTime << " s\n";
        }
}

// -------------------------------------
// Reader thread for k > 10
// -------------------------------------
struct ReaderArgs {
    std::string filename;
    std::streampos startPos;
    std::streampos endPos;
    int k;
    std::vector<PackageManager>* pms;
    int numPMs;
    std::atomic<uint64_t>* totalKmers;
    int readerId;
    uint64_t totalLengthProcessed = 0;
    uint64_t totalSequencesProcessed = 0;
};

static void readerFunc(ReaderArgs args)
{
    auto start = std::chrono::steady_clock::now();

    double currentBlockScore = 0.0;
    std::ifstream mafFile(args.filename);
    mafFile.seekg(args.startPos);
    uint64_t localCount = 0;

    std::vector<BinBuffer> localBins(1 << P1_BITS);
    for (size_t i = 0; i < localBins.size(); ++i)
        localBins[i].p1_id = static_cast<uint16_t>(i);

    /*  mask keeps only the low 2*k bits after we shift‑in a new base        *
     *  (If k == 64 the whole 128‑bit word is used, so we just take ~0ULL).  */
    const uint128_t MASK = (args.k == 64)
        ? std::numeric_limits<uint128_t>::max()
        : (((uint128_t)1 << (2 * args.k)) - 1);

    std::string line;
    while (true)
    {
        if (mafFile.tellg() >= args.endPos) break;
        if (!std::getline(mafFile, line))   break;
        if (line.empty())                   continue;

        // Block‑header line ― record its score.
        if (line[0] == 'a')
        {
            std::size_t pos = line.find("score=");
            if (pos != std::string::npos)
                try { currentBlockScore = std::stod(line.substr(pos + 6)); }
                catch (...) { currentBlockScore = 0.0; }
            else
                currentBlockScore = 0.0;
            continue;
        }

        if (line[0] != 's') continue;               // we only want "s" lines
        if (currentBlockScore < g_minAScore ||
            currentBlockScore > g_maxAScore)        // A‑score filter
            continue;

        // -------- parse the "s" line --------
        std::istringstream iss(line);
        std::string token, src, seq;
        iss >> token >> src;
        for (int i = 0; i < 4; ++i) iss >> token;   // start, size, strand, srcSize
        iss >> seq;                                 // aligned sequence (with dashes)

        // Genome name → numeric id
        size_t dotPos = src.find('.');
        std::string gname = (dotPos == std::string::npos)
                            ? src
                            : src.substr(0, dotPos);
        if (g_useGenomeIdFilter &&
            g_includedGenomeNames.find(gname) == g_includedGenomeNames.end())
            continue;
        uint8_t gID = getGenomeID(gname);

        // Optional quality line
        std::string cleanedQ;
        if ((g_minQVal >= 0 || g_maxQVal >= 0))
        {
            std::streampos savePos = mafFile.tellg();
            std::string qline;
            if (std::getline(mafFile, qline))
            {
                if (!qline.empty() && qline[0] == 'q')
                {
                    std::istringstream qiss(qline);
                    std::string qtoken, qsrc, qseq;
                    qiss >> qtoken >> qsrc;
                    for (int i = 0; i < 4; ++i) qiss >> qtoken;
                    qiss >> qseq;
                    cleanedQ = removeDashes(qseq);      // still using existing helper
                    std::transform(cleanedQ.begin(), cleanedQ.end(), cleanedQ.begin(),
               [](char c){ return std::toupper(c); });
                }
                else
                {
                    mafFile.seekg(savePos);
                }
            }
        }

        // Remove gaps from the sequence  (unchanged helper)
        std::string cleaned = removeDashes(seq);
        std::transform(cleaned.begin(), cleaned.end(), cleaned.begin(),
        [](char c){ return std::toupper(c); });
        args.totalLengthProcessed    += cleaned.length();
        args.totalSequencesProcessed += 1;

        if (cleaned.size() < static_cast<size_t>(args.k)) continue;

        /* -----------------------------------------------------------------
         *  ROLLING  K‑MER  LOOP
         * ----------------------------------------------------------------*/
        uint128_t kmerBits = 0;      // accumulates the forward strand
        int       validLen = 0;      // # of consecutive valid (non‑N, pass‑Q) bases

        for (size_t i = 0; i < cleaned.size(); ++i)
        {
            char c = cleaned[i];
            uint8_t val;
            switch (c)               // convert DNA to 2 bits
            {
                case 'A': case 'a': val = 0; break;
                case 'C': case 'c': val = 1; break;
                case 'G': case 'g': val = 2; break;
                case 'T': case 't': val = 3; break;
                default:            // any non‑ACGT breaks the window
                    validLen = 0;
                    kmerBits = 0;
                    continue;
            }

            // Quality check (if requested)
            if (!cleanedQ.empty() && !isQValid(cleanedQ[i]))
            {
                validLen = 0;
                kmerBits = 0;
                continue;
            }

            /* shift the previous k‑mer left by 2 bits and add the
             * new base; then mask to keep at most 2*k bits           */
            kmerBits = ((kmerBits << 2) | val) & MASK;

            if (++validLen < args.k) continue;       // not enough bases yet

            uint128_t kmerOut = kmerBits;
            if (g_useCanonical)
                kmerOut = canonicalKmer128(kmerOut, args.k);   // unchanged helper

            uint16_t p1 = extractP1(kmerOut, args.k);
            auto &vec   = localBins[p1].records;
            vec.push_back({kmerOut, gID});
            ++localCount;

            if (vec.size() >= BIN_CAPACITY)          // flush full bin
            {
                BinChunk chunk = partialSortAndCompact(localBins[p1], args.k);
                vec.clear();
                int pmIndex = static_cast<int>(chunk.p1_id % args.numPMs);
                auto &pm = (*args.pms)[pmIndex];
                {
                    std::lock_guard<std::mutex> lk(pm.pmMutex);
                    pm.pmQueue.push(std::move(chunk));
                    pm.pmCV.notify_one();
                }
            }
        }
    }

    /* flush any leftover records in the per‑thread bins */
    for (auto &b : localBins)
        if (!b.records.empty())
        {
            BinChunk chunk = partialSortAndCompact(b, args.k);
            b.records.clear();
            int pmIndex = static_cast<int>(chunk.p1_id % args.numPMs);
            auto &pm = (*args.pms)[pmIndex];
            {
                std::lock_guard<std::mutex> lk(pm.pmMutex);
                pm.pmQueue.push(std::move(chunk));
                pm.pmCV.notify_one();
            }
        }

    args.totalKmers->fetch_add(localCount, std::memory_order_relaxed);

    auto end = std::chrono::steady_clock::now();
    double duration = std::chrono::duration<double>(end - start).count();
    {
        std::lock_guard<std::mutex> lock(coutMutex);
        std::cout << "Reader thread " << args.readerId
                  << " finished in " << duration << " seconds.\n";
    }

    if (args.totalSequencesProcessed > 0)
    {
        std::cout << "Avg sequence size processed: "
                  << args.totalLengthProcessed / args.totalSequencesProcessed
                  << "\n";
    }
}

// ---------------------------------------------
// Merge multiple bin_X_*.bin => bin_X_complete => final.bin
// ---------------------------------------------
struct FileHeader {
    char     magic[4];
    uint32_t k;
    uint16_t p1_bits, p2_bits;
    uint8_t  genome_id_bits, suffix_mode, reserved[2];
};

struct GroupMeta { uint16_t p2_val; uint32_t sz; };

/* Combined container */
struct Comb {
    uint8_t mode{};
    std::unordered_map<uint16_t,std::vector<uint32_t>> v0;
    std::unordered_map<uint16_t,std::vector<uint64_t>> v1;
    std::unordered_map<uint16_t,std::vector<uint128_t>> v3;
};


// Remove all "bin_" files if we want total purge
static void purgeIntermediateFiles(const std::string &excludeFile = "")
{
    for (auto &entry : std::filesystem::directory_iterator(g_tempFilesDir)) {
        if (!entry.is_regular_file())
            continue;
        std::string fname = entry.path().filename().string();
        if (fname.rfind("bin_", 0) == 0) {
            if (fname == excludeFile)
                continue;
            std::error_code ec;
            std::filesystem::remove(entry.path(), ec);
            if (ec) {
                std::cerr << "Warning: could not remove " << fname << ": " << ec.message() << "\n";
            }
        }
    }
}



static std::vector<char> buildPayload(
    uint16_t p1,
    const std::vector<std::filesystem::path>& files,
    int mode,
    uint32_t &k_out,
    double &readTime,
    double &sortTime
) {
    Comb comb; comb.mode = (uint8_t)mode;
    uint32_t k_val = 0;
    bool first = true;

    auto t0 = std::chrono::steady_clock::now();
    // read all bin_<p1>_*.bin fragments
    for (auto &path : files) {
        std::ifstream ifs(path, std::ios::binary);
        if (!ifs) continue;
        FileHeader hdr;
        ifs.read(reinterpret_cast<char*>(&hdr), sizeof(hdr));
        if (std::memcmp(hdr.magic, "KMHD", 4) != 0 ||
            hdr.suffix_mode != comb.mode)
            continue;
        if (first) { k_val = hdr.k; first = false; }
        // skip sentinel
        ifs.seekg(4, std::ios::cur);
        uint32_t gc; ifs.read(reinterpret_cast<char*>(&gc), 4);
        std::vector<GroupMeta> gmeta(gc);
        for (auto &gm : gmeta) {
            ifs.read(reinterpret_cast<char*>(&gm.p2_val), 2);
            ifs.read(reinterpret_cast<char*>(&gm.sz),     4);
        }
        // skip sentinel
        ifs.seekg(4, std::ios::cur);

        // read each group
        for (auto &gm : gmeta) {
            if (!gm.sz) continue;
            switch (mode) {
            case 0: {
                auto &v = comb.v0[gm.p2_val];
                size_t old = v.size();
                v.resize(old + gm.sz);
                ifs.read(reinterpret_cast<char*>(v.data()+old),
                         gm.sz * sizeof(uint32_t));
            } break;
            case 1: {
                auto &v = comb.v1[gm.p2_val];
                size_t old = v.size();
                v.resize(old + gm.sz);
                ifs.read(reinterpret_cast<char*>(v.data()+old),
                         gm.sz * sizeof(uint64_t));
            } break;
            case 3: {
                auto &v = comb.v3[gm.p2_val];
                size_t old = v.size();
                v.resize(old + gm.sz);
                for (size_t i = 0; i < gm.sz; ++i) {
                    unsigned char buf[16];
                    ifs.read(reinterpret_cast<char*>(buf), 16);
                    uint128_t x = 0;
                    for (int b = 0; b < 16; ++b)
                        x |= (uint128_t)buf[b] << (8*b);
                    v[old + i] = x;
                }
            } break;
            }
        }
    }
    auto t1 = std::chrono::steady_clock::now();
    readTime += std::chrono::duration<double>(t1 - t0).count();

    // sort
    auto t2 = std::chrono::steady_clock::now();
    if (mode == 0) for (auto &kv : comb.v0) LSD32(kv.second);
    if (mode == 1) for (auto &kv : comb.v1) LSD64(kv.second);
    if (mode == 3) for (auto &kv : comb.v3) LSD128(kv.second);
    auto t3 = std::chrono::steady_clock::now();
    sortTime += std::chrono::duration<double>(t3 - t2).count();

    // collect sorted p2 keys
    std::vector<uint16_t> p2vals;
    if (mode == 0) for (auto &kv : comb.v0) p2vals.push_back(kv.first);
    if (mode == 1) for (auto &kv : comb.v1) p2vals.push_back(kv.first);
    if (mode == 3) for (auto &kv : comb.v3) p2vals.push_back(kv.first);
    std::sort(p2vals.begin(), p2vals.end());

    // serialize into a flat buffer
    std::vector<char> out;
    auto W = [&](const void *ptr, size_t n){
        const char *c = reinterpret_cast<const char*>(ptr);
        out.insert(out.end(), c, c + n);
    };

    FileHeader hdr;
    std::memcpy(hdr.magic, "KMHD", 4);
    hdr.k            = k_val;
    hdr.p1_bits      = P1_BITS;
    hdr.p2_bits      = P2_BITS;
    hdr.genome_id_bits = 8;
    hdr.suffix_mode  = (uint8_t)mode;
    hdr.reserved[0] = hdr.reserved[1] = 0;
    W(&hdr, sizeof(hdr));
    W("kmhd", 4);

    uint32_t gc = p2vals.size();
    W(&gc, 4);
    for (auto p2 : p2vals) {
        W(&p2, 2);
        uint32_t sz =
          (mode==0 ? comb.v0[p2].size()
         : mode==1 ? comb.v1[p2].size()
                   : comb.v3[p2].size());
        W(&sz, 4);
    }
    W("kmhd", 4);

    for (auto p2 : p2vals) {
        if (mode == 0) {
            auto &v = comb.v0[p2];
            W(v.data(), v.size()*4);
        } else if (mode == 1) {
            auto &v = comb.v1[p2];
            W(v.data(), v.size()*8);
        } else {
            for (auto &x : comb.v3[p2]) {
                unsigned char buf[16];
                for (int j = 0; j < 16; ++j)
                    buf[j] = (unsigned char)((x >> (8*j)) & 0xFF);
                W(buf,16);
            }
        }
    }

    k_out = k_val;
    return out;
}

// Main
int main(int argc, char* argv[])
{
    // These will hold the user-provided values for thread parameters
    int userProvidedReaders = -1;
    int userProvidedPMs = -1;
    int userProvidedThreads = -1;

    std::vector<std::string> positionalArgs;

    for (int i = 1; i < argc; i++) {
        std::string arg(argv[i]);
        if (arg == "-c") {
            g_useCanonical = true;
        }
        else if (arg.rfind("--purge_intermediate", 0) == 0) {
            size_t pos = arg.find('=');
            if (pos != std::string::npos) {
                std::string value = arg.substr(pos + 1);
                g_purgeIntermediate = (value == "true" || value == "1");
            } else {
                g_purgeIntermediate = true;
            }
        }
        else if (arg == "--binary_file_output") {
            if (i + 1 < argc) {
                g_binaryOutputFile = argv[++i];
            } else {
                std::cerr << "Missing argument for --binary_file_output\n";
                return 1;
            }
        }
        else if (arg.rfind("--genome_ids", 0) == 0) {
            size_t pos = arg.find('=');
            if (pos != std::string::npos) {
                std::string val = arg.substr(pos + 1);
                if (val != "all") {
                    g_useGenomeIdFilter = true;
                    std::stringstream ss(val);
                    std::string token;
                    while (std::getline(ss, token, ',')) {
                        if (!token.empty()) {
                            g_includedGenomeNames.insert(token);
                        }
                    }
                }
            }
        }
        else if (arg.rfind("--min_a_score", 0) == 0) {
            size_t pos = arg.find('=');
            if (pos != std::string::npos) {
                std::string val = arg.substr(pos + 1);
                g_minAScore = std::stod(val);
            }
        }
        else if (arg.rfind("--max_a_score", 0) == 0) {
            size_t pos = arg.find('=');
            if (pos != std::string::npos) {
                std::string val = arg.substr(pos + 1);
                g_maxAScore = std::stod(val);
            }
        }
        else if (arg.rfind("--min_q_level", 0) == 0) {
            size_t pos = arg.find('=');
            if (pos != std::string::npos) {
                std::string val = arg.substr(pos + 1);
                if (!val.empty()) {
                    int parsed = -1;
                    if (val.size() == 1) {
                        parsed = parseQVal(val[0]);
                    }
                    if (parsed >= 0) {
                        g_minQVal = parsed;
                    } else {
                        std::cerr << "Invalid value for --min_q_level\n";
                        return 1;
                    }
                }
            }
        }
        else if (arg.rfind("--max_q_level", 0) == 0) {
            size_t pos = arg.find('=');
            if (pos != std::string::npos) {
                std::string val = arg.substr(pos + 1);
                if (!val.empty()) {
                    int parsed = -1;
                    if (val.size() == 1) {
                        parsed = parseQVal(val[0]);
                    }
                    if (parsed >= 0) {
                        g_maxQVal = parsed;
                    } else {
                        std::cerr << "Invalid value for --max_q_level\n";
                        return 1;
                    }
                }
            }
        }
        else if (arg == "--k") {
            if (i + 1 < argc) {
                g_k = std::stoi(argv[++i]);
            } else {
                std::cerr << "Missing argument for --k\n";
                return 1;
            }
        }
        else if (arg == "--reader_threads") {
            if (i + 1 < argc) {
                userProvidedReaders = std::stoi(argv[++i]);
            } else {
                std::cerr << "Missing argument for --reader_threads\n";
                return 1;
            }
        }
        else if (arg == "--package_manager_threads") {
            if (i + 1 < argc) {
                userProvidedPMs = std::stoi(argv[++i]);
            } else {
                std::cerr << "Missing argument for --package_manager_threads\n";
                return 1;
            }
        }
        else if (arg == "--threads") {
            if (i + 1 < argc) {
                userProvidedThreads = std::stoi(argv[++i]);
            } else {
                std::cerr << "Missing argument for --threads\n";
                return 1;
            }
        }
        else if (arg == "--temp_files_dir") {
            if (i + 1 < argc) {
                g_tempFilesDir = argv[++i];
            } else {
                std::cerr << "Missing argument for --temp_files_dir\n";
                return 1;
            }
        }
        else if (arg == "--output_directory") {
            if (i + 1 < argc) {
                g_outputDirectory = argv[++i];
            } else {
                std::cerr << "Missing argument for --output_directory\n";
                return 1;
            }
        }
        else {
            positionalArgs.push_back(arg);
        }
    }

    if (!std::filesystem::exists(g_tempFilesDir)) {
        if (!std::filesystem::create_directories(g_tempFilesDir)) {
            std::cerr << "Error: unable to create temporary files directory: " << g_tempFilesDir << "\n";
            return 1;
        }
    }
    if (!std::filesystem::exists(g_outputDirectory)) {
        if (!std::filesystem::create_directories(g_outputDirectory)) {
            std::cerr << "Error: unable to create output directory: " << g_outputDirectory << "\n";
            return 1;
        }
    }

    if (positionalArgs.empty()) {
        std::cerr << "Usage:\n"
                  << "  " << argv[0] << " [ -c ]\n"
                  << "        [--purge_intermediate[=true|false]]\n"
                  << "        [--binary_file_output <filename>]\n"
                  << "        [--genome_ids=all|name1,name2,...]\n"
                  << "        [--min_a_score=<val>] [--max_a_score=<val>]\n"
                  << "        [--min_q_level=0..9|F] [--max_q_level=0..9|F]\n"
                  << "        --k <VAL>\n"
                  << "        [--reader_threads <VAL>] [--package_manager_threads <VAL>] | [--threads <VAL>]\n"
                  << "        [--temp_files_dir <DIR>] [--output_directory <DIR>]\n"
                  << "        <MAFfile>\n\n"
                  << "Notes:\n"
                  << "  --k is required. Threads can be specified either by:\n"
                  << "       --reader_threads and --package_manager_threads (both)\n"
                  << "       OR\n"
                  << "       --threads (which is then split ~2:1 for readers:PMs).\n"
                  << "  --temp_files_dir is where intermediate bin files go (default: current dir).\n"
                  << "  --output_directory is where final.bin and final.metadata go (default: current dir).\n"
                  << std::endl;
        return 1;
    }

    std::string mafFileOriginal = positionalArgs[0];

    if (userProvidedThreads != -1 && (userProvidedReaders != -1 || userProvidedPMs != -1)) {
        std::cerr << "Error: --threads is mutually exclusive with --reader_threads and --package_manager_threads\n";
        return 1;
    }
    else if (userProvidedThreads != -1) {
        if (userProvidedThreads <= 0) {
            std::cerr << "Error: --threads must be > 0\n";
            return 1;
        }
        g_numReaders = 0.3 * userProvidedThreads; 
        if (g_numReaders < 1) g_numReaders = 1;
        g_numPMs = userProvidedThreads - g_numReaders; 
        if (g_numPMs < 1) g_numPMs = 1;
    }
    else {
        if (userProvidedReaders <= 0 || userProvidedPMs <= 0) {
            std::cerr << "Error: must provide --reader_threads <val> and --package_manager_threads <val>, or use --threads\n";
            return 1;
        }
        g_numReaders = userProvidedReaders;
        g_numPMs = userProvidedPMs;
    }

    if (g_k <= 0) {
        std::cerr << "Error: k not specified (or invalid). Use --k <VAL>\n";
        return 1;
    }
    if (g_k > KMER_MAX_LENGTH) {
        std::cerr << "k too large (max " << KMER_MAX_LENGTH << ")\n";
        return 1;
    }

    bool isGZ = endsWithGZ(mafFileOriginal);
    std::string decompressedFile;
    std::string mafFileToProcess = mafFileOriginal;
    if (isGZ) {
        decompressedFile = decompressGZ(mafFileOriginal);
        if (decompressedFile.empty()) {
            std::cerr << "Decompression failed. Exiting.\n";
            return 1;
        }
        mafFileToProcess = decompressedFile;
    }

    if (g_purgeIntermediate) {
        purgeIntermediateFiles(g_binaryOutputFile);
    }

    std::vector<std::streampos> chunkStarts;
    findMAFChunkBoundaries(mafFileToProcess, g_numReaders, chunkStarts);
    std::atomic<uint64_t> totalKmers(0);

    // If k <= 10 => smallK approach
    if (g_k <= 10) {
        std::vector<std::thread> readers;
        std::vector<LocalHashData> localData(g_numReaders);

        for (int i = 0; i < g_numReaders; i++) {
            SmallKReaderArgs rargs {
                mafFileToProcess,
                chunkStarts[i],
                chunkStarts[i+1],
                g_k,
                &localData[i],
                &totalKmers,
                i
            };
            readers.emplace_back(smallKReaderFunc, rargs);
        }
        for (auto &t : readers) {
            t.join();
        }

        std::cout << "Total kmers processed: " << totalKmers.load() << "\n";

        google::dense_hash_map<uint32_t, std::array<uint32_t, 256>> globalMap;
        globalMap.set_empty_key(0xFFFFFFFF);
        globalMap.set_deleted_key(0xFFFFFFFE);

        for (int i = 0; i < g_numReaders; i++) {
            for (auto &kv : localData[i].map) {
                auto it = globalMap.find(kv.first);
                if (it == globalMap.end()) {
                    globalMap.insert({kv.first, kv.second});
                } else {
                    for (int g = 0; g < 256; g++) {
                        it->second[g] += kv.second[g];
                    }
                }
            }
        }

        auto finalStart = std::chrono::steady_clock::now();
        std::string finalBinPath = g_outputDirectory + "/" + g_binaryOutputFile;
        std::ofstream ofs(finalBinPath, std::ios::binary);
        if (!ofs) {
            std::cerr << "Error creating " << finalBinPath << "\n";
            if (isGZ) removeFileIfExists(decompressedFile);
            return 0;
        }
        ofs.write("BINF", 4);
        ofs.write((const char*)&g_k, sizeof(g_k));
        uint8_t nIDs = (uint8_t)g_genomeNames.size();
        ofs.write((const char*)&nIDs, sizeof(nIDs));
        for (uint8_t i = 0; i < nIDs; i++) {
            const std::string &nm = g_genomeNames[i];
            uint16_t len = (uint16_t)nm.size();
            ofs.write((const char*)&len, sizeof(len));
            ofs.write(nm.data(), len);
        }
        ofs.write("ENDG", 4);
        uint64_t mapSize = (uint64_t)globalMap.size();
        ofs.write((const char*)&mapSize, sizeof(mapSize));
        for (auto &kv : globalMap) {
            uint32_t key32 = kv.first;
            ofs.write((const char*)&key32, sizeof(key32));
            ofs.write((const char*)kv.second.data(), 256 * sizeof(uint32_t));
        }
        ofs.write("ENDM", 4);
        ofs.close();

        auto finalEnd = std::chrono::steady_clock::now();
        double finalDuration = std::chrono::duration<double>(finalEnd - finalStart).count();
        {
            std::lock_guard<std::mutex> lock(coutMutex);
            std::cout << "Final.bin creation took " << finalDuration << " seconds.\n";
        }
        std::cout << finalBinPath << " created.\n";
    }
    else {
        // For k > 10
        std::vector<PackageManager> pms(g_numPMs);
        std::vector<std::thread> pmThreads;
        for (int i = 0; i < g_numPMs; i++) {
            pmThreads.emplace_back(PMThreadLoop, std::ref(pms[i]), g_k, i);
        }
        std::vector<std::thread> readers;
        for (int i = 0; i < g_numReaders; i++) {
            ReaderArgs rargs {
                mafFileToProcess, chunkStarts[i], chunkStarts[i+1],
                g_k, &pms, g_numPMs, &totalKmers,
                i
            };
            readers.emplace_back(readerFunc, rargs);
        }
        for (auto &t : readers) {
            t.join();
        }
        for (auto &pm : pms) {
            std::lock_guard<std::mutex> lk(pm.pmMutex);
            pm.pmReadersDone = true;
            pm.pmCV.notify_one();
        }
        for (auto &t : pmThreads) {
            t.join();
        }
        std::cout << "Aggregating done. Total kmers: " << totalKmers.load() << "\n";

        
        size_t maxId = 0;
        for (auto &kv : g_genomeMap) maxId = std::max(maxId, (size_t)kv.second);
        std::vector<std::string> gNames(maxId + 1);
        for (auto &kv : g_genomeMap) gNames[kv.second] = kv.first;

        

        std::vector<uint16_t> allBins;
        allBins.reserve(1 << P1_BITS);
        for (uint16_t i = 0; i < (1 << P1_BITS); i++) {
            allBins.push_back(i);
        }
        std::vector<std::string> allFiles;
        for (auto &entry : std::filesystem::directory_iterator(g_tempFilesDir)) {
            if (!entry.is_regular_file()) continue;
            std::string fname = entry.path().filename().string();
            if (fname.rfind("bin_", 0) == 0) {
                allFiles.push_back(entry.path().string());
            }
        }
        std::unordered_map<uint16_t, std::vector<std::string>> binMap;
        auto parseP1 = [&](const std::string &fn) -> uint16_t {
            // We expect something like ".../bin_<p1>_<rnd>.bin"
            std::filesystem::path p(fn);
            std::string base = p.filename().string(); 
            size_t pos = base.find('_');
            if (pos == std::string::npos) return 0xFFFF;
            size_t pos2 = base.find('_', pos + 1);
            if (pos2 == std::string::npos) return 0xFFFF;
            std::string p1str = base.substr(pos + 1, pos2 - (pos + 1));
            return (uint16_t)std::stoi(p1str);
        };
        for (auto &f : allFiles) {
            uint16_t p1val = parseP1(f);
            if (p1val != 0xFFFF) {
                binMap[p1val].push_back(f);
            }
        }

        std::mutex outMutex;
        size_t nextBinIdx = 0;
        std::vector<std::thread> binThreads;
        size_t threadCount = (g_numPMs + g_numReaders);
        if (threadCount < 1) threadCount = 1;

      
        std::unordered_map<uint16_t, std::vector<std::filesystem::path>> bins;
        {
            std::regex re(R"(bin_(\d+)_.*\.bin)");
            for (auto &e : std::filesystem::directory_iterator(g_tempFilesDir)) {
                if (!e.is_regular_file()) continue;
                std::smatch m;
                auto fn = e.path().filename().string();
                if (std::regex_match(fn, m, re)) {
                    uint16_t p1 = static_cast<uint16_t>(std::stoi(m[1].str()));
                    bins[p1].push_back(e.path());
                }
            }
        }

        // 3) Create and write the global header of final.bin
        const std::string finalBinPath = g_outputDirectory + "/" + g_binaryOutputFile;
        {
            int fd = ::open(finalBinPath.c_str(),
                            O_RDWR | O_CREAT | O_TRUNC,
                            0666);
            if (fd < 0) {
                std::cerr << "Error: cannot create " << finalBinPath << "\n";
                return 1;
            }

            // write "BINF"
            size_t offs = 0;
            auto PW = [&](const void* p, size_t n) {
                pwrite(fd, p, n, offs);
                offs += n;
            };

            PW("BINF", 4);
            // placeholder for k; we’ll patch this once
            uint32_t k_placeholder = 0;
            PW(&k_placeholder, sizeof(k_placeholder));

            // genome‐name table
            uint8_t nIDs = (uint8_t)gNames.size();
            PW(&nIDs, 1);
            for (auto &nm : gNames) {
                uint16_t L = (uint16_t)nm.size();
                PW(&L, 2);
                PW(nm.data(), L);
            }
            PW("ENDG", 4);
            ::close(fd);
        }

        // track where the next packet can go
        struct stat st; stat(finalBinPath.c_str(), &st);
        std::atomic<uint64_t> nextOff(st.st_size), fileSz(st.st_size);
        std::mutex offMu, metaMu;
        std::ofstream metaOut(g_outputDirectory + "/final.metadata",
                              std::ios::trunc);

        // make a sorted list of p1 IDs to process
        std::vector<uint16_t> p1list;
        p1list.reserve(bins.size());
        for (auto &kv : bins) p1list.push_back(kv.first);
        std::sort(p1list.begin(), p1list.end());

        std::vector<uint16_t> normals, outliers;
        for (auto p1 : p1list) {
            if (outlierBins.count(p1)) outliers.push_back(p1);
            else                     normals.push_back(p1);
        }

        // interleave: after every N normals, insert one outlier
        std::vector<uint16_t> interleaved;
        size_t normalPerOut = normals.size() / (outliers.empty() ? 1 : outliers.size());
        size_t ni = 0, oi = 0;
        while (ni < normals.size() || oi < outliers.size()) {
            for (size_t k = 0; k < normalPerOut && ni < normals.size(); ++k)
                interleaved.push_back(normals[ni++]);
            if (oi < outliers.size())
                interleaved.push_back(outliers[oi++]);
        }
        while (ni < normals.size()) interleaved.push_back(normals[ni++]);

        p1list.swap(interleaved);

        std::atomic<size_t> idx(0);
        if (threadCount == 0) threadCount = 1;

        // 4) Worker that for each p1:
        //    a) merges its bin files, sorts them via buildPayload()
        //    b) patches global k into header once
        //    c) reserves space in final.bin, writes "BINc", p1, sendfile payload, "BEND"
        //    d) records offset in final.metadata, and asynchronously deletes temp files
        std::vector<TimeStats> threadStats(threadCount);
        auto worker = [&](int tid){
            double read_time = 0, sort_time = 0, write_time = 0;
            bool patchedK = false;
        
            int remainder_bits = 2 * g_k - P1_BITS - P2_BITS;
            int payloadMode   = determineSuffixMode(remainder_bits);
        
            while (true) {
                size_t i = idx.fetch_add(1);
                if (i >= p1list.size()) break;
                uint16_t p1 = p1list[i];
                auto &files = bins[p1];
                uint32_t local_k = 0;
        
                bool isOut = outlierBins.count(p1) > 0;
                // if this is a slow bin, grab the mutex so no other thread
                // can sort/write another outlier at the same time
                if (isOut) {
                    std::lock_guard<std::mutex> L(outlierMutex);
                    // fall through into the normal processing code
                }
        
                // a) build & sort payload
                auto payload = buildPayload(p1, files, payloadMode,
                                            local_k, read_time, sort_time);
        
                // b) patch global k (only once)
                if (!patchedK) {
                    int fd = ::open(finalBinPath.c_str(), O_WRONLY);
                    pwrite(fd, &local_k, sizeof(local_k), 4);
                    ::close(fd);
                    patchedK = true;
                }
        
                // c) reserve space
                uint64_t segSz = 4 + 2 + payload.size() + 4;
                uint64_t myOff;
                {
                    std::lock_guard<std::mutex> lk(offMu);
                    myOff = nextOff.fetch_add(segSz);
                    uint64_t cur = nextOff.load();
                    if (cur > fileSz.load()) {
                        fileSz.store(cur);
                        int ffd = ::open(finalBinPath.c_str(), O_WRONLY);
                        ftruncate(ffd, cur);
                        ::close(ffd);
                    }
                }
                // record in metadata
                {
                    std::lock_guard<std::mutex> lk(metaMu);
                    metaOut << p1 << " " << myOff << "\n";
                }
        
                // d) write payload -> final.bin using std::ofstream
                {
                    auto t0 = std::chrono::steady_clock::now();
                    std::ofstream ofs(finalBinPath,
                                      std::ios::in | std::ios::out | std::ios::binary);
                    if (!ofs) {
                        std::cerr << "Error opening " << finalBinPath << " for writing\n";
                        std::exit(1);
                    }
                    ofs.seekp(myOff);
                    ofs.write("BINc", 4);
                    ofs.write(reinterpret_cast<const char*>(&p1), sizeof(p1));
                    ofs.write(payload.data(), payload.size());
                    ofs.write("BEND", 4);
                    ofs.flush();
                    auto t1 = std::chrono::steady_clock::now();
                    write_time += std::chrono::duration<double>(t1 - t0).count();
                }
        
                // e) asynchronously delete the just-written bin files
                std::vector<std::filesystem::path> to_remove = files;
                std::thread([to_remove]() {
                    for (auto &path : to_remove) {
                        std::error_code ec;
                        std::filesystem::remove(path, ec);
                    }
                }).detach();
            }
        
            threadStats[tid] = TimeStats{read_time, sort_time, write_time};
        };
        

        // 5) Launch writers
        std::vector<std::thread> thr;
        for (unsigned t = 0; t < threadCount; ++t)
            thr.emplace_back(worker, t);
        for (auto &t : thr) t.join();
        metaOut.close();

        // 6) (Optional) print per-thread timings
        for (unsigned t = 0; t < threadCount; ++t) {
            auto &s = threadStats[t];
            std::cout << "[TIMING][T" << t << "] "
                      << "read="  << s.readTime  << "s, "
                      << "sort="  << s.sortTime  << "s, "
                      << "write=" << s.writeTime << "s\n";
        }

        std::cout << "Done. final.bin written to " << finalBinPath << "\n";
    }

    return 0;
}
