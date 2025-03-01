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
#include <google/dense_hash_map>
#include <unordered_set>
#include <limits>
#include <zlib.h>
#include <fcntl.h>      // for open
#include <sys/mman.h>   // for mmap, munmap
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

// ---------------------- NEW: decodeKmer for 128 bits ----------------------
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

// ---------------------- NEW: Reverse complement logic (128-bit) ----------------------
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
                } else {
                    mafFile.seekg(savePos);
                }
            }
        }

        std::string cleaned = removeDashes(seq);
        
        if (cleaned.size() < (size_t)args.k) continue;

        for (size_t pos = 0; pos + args.k <= cleaned.size(); pos++) {
            uint32_t kmer32 = 0;
            bool valid = true;
            for (int i = 0; i < args.k; i++) {
                char c = cleaned[pos + i];
                uint32_t val;
                switch (c) {
                    case 'A': val = 0; break;
                    case 'C': val = 1; break;
                    case 'G': val = 2; break;
                    case 'T': val = 3; break;
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
// K > 10 aggregator approach (128-bit kmers, extended to handle up to 70)
// ---------------------------------------------
struct KmerRecord {
    uint128_t kmer;
    uint8_t   genomeID;
};

// For remainder bits extraction
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
    // 0 => up to 32 bits suffix
    // 1 => up to 64 bits suffix
    // 3 => up to 120 bits suffix + 8 bits genome ID (total 128)
    uint8_t suffixMode;
    std::vector<uint32_t> mode0;    // remainder<<8 | genomeID
    std::vector<uint64_t> mode1;    // remainder<<8 | genomeID
    std::vector<uint128_t> mode3;   // remainder<<8 | genomeID
};

struct BinChunk {
    uint16_t p1_id;
    std::vector<SuffixGroup> groups;
};

static void LSDRadixSort32(std::vector<uint32_t> &arr)
{
    const int PASSES = 4;
    const int RADIX = 256;
    std::vector<uint32_t> tmp(arr.size());
    for (int pass = 0; pass < PASSES; pass++){
        int shift = pass * 8;
        uint32_t count[RADIX];
        std::memset(count, 0, sizeof(count));
        for (auto v : arr) {
            uint8_t b = (uint8_t)((v >> shift) & 0xFF);
            count[b]++;
        }
        uint32_t pos[RADIX];
        pos[0] = 0;
        for (int i = 1; i < RADIX; i++) {
            pos[i] = pos[i - 1] + count[i - 1];
        }
        for (auto v : arr) {
            uint8_t b = (uint8_t)((v >> shift) & 0xFF);
            tmp[pos[b]++] = v;
        }
        arr.swap(tmp);
    }
}

static void LSDRadixSort64(std::vector<uint64_t> &arr)
{
    const int PASSES = 8;
    const int RADIX = 256;
    std::vector<uint64_t> tmp(arr.size());
    for (int pass = 0; pass < PASSES; pass++){
        int shift = pass * 8;
        uint64_t count[RADIX];
        std::memset(count, 0, sizeof(count));
        for (auto v : arr) {
            uint8_t b = (uint8_t)((v >> shift) & 0xFF);
            count[b]++;
        }
        uint64_t pos[RADIX];
        pos[0] = 0;
        for (int i = 1; i < RADIX; i++) {
            pos[i] = pos[i - 1] + count[i - 1];
        }
        for (auto v : arr) {
            uint8_t b = (uint8_t)((v >> shift) & 0xFF);
            tmp[pos[b]++] = v;
        }
        arr.swap(tmp);
    }
}

// NEW: LSDRadixSort for 128-bit values
static void LSDRadixSort128(std::vector<uint128_t> &arr)
{
    const int PASSES = 16; 
    const int RADIX = 256;
    std::vector<uint128_t> tmp(arr.size());
    for (int pass = 0; pass < PASSES; pass++){
        int shift = pass * 8;
        uint64_t count[RADIX];
        std::memset(count, 0, sizeof(count));
        for (auto &v : arr) {
            // Extracting the byte:
            uint128_t b = (v >> shift) & 0xFFU;
            count[b.convert_to<uint64_t>()]++;
        }
        uint64_t pos[RADIX];
        pos[0] = 0;
        for (int i = 1; i < RADIX; i++) {
            pos[i] = pos[i - 1] + count[i - 1];
        }
        for (auto &v : arr) {
            uint128_t b = (v >> shift) & 0xFFU;
            uint64_t idx = b.convert_to<uint64_t>();
            tmp[pos[idx]++] = v;
        }
        arr.swap(tmp);
    }
}

// Decide suffix mode based on remainder bits
//  remainder_bits <= 32 => mode 0
//  remainder_bits <= 64 => mode 1
//  remainder_bits <=120 => mode 3
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
    // All temp (intermediate) files now go into g_tempFilesDir
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

    // In this aggregated file, all groups have the same suffix mode
    // but let's just pick the first group's mode if it exists
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
                // Write each 128-bit value in little-endian
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
                writePackageToDisk(c.p1_id, acc, k);
            }
            lock.lock();
        }
        if (pm.pmReadersDone && pm.pmQueue.empty()) break;
    }
    for (auto &kv : pm.pmData) {
        if (kv.second.totalCount > 0) {
            writePackageToDisk(kv.first, kv.second, k);
        }
    }
    auto end = std::chrono::steady_clock::now();
    double duration = std::chrono::duration<double>(end - start).count();
    {
        std::lock_guard<std::mutex> lock(coutMutex);
        std::cout << "PM thread " << pmId << " finished in " << duration << " seconds.\n";
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
    for (size_t i = 0; i < localBins.size(); i++) {
        localBins[i].p1_id = (uint16_t)i;
    }

    std::string line;
    while (true) {
        if (mafFile.tellg() >= args.endPos) break;
        if (!std::getline(mafFile, line)) break;
        if (line.empty()) continue;

        if (line[0] == 'a') {
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
        for(int i = 0; i < 4; i++) iss >> token;
        iss >> seq;

        size_t dotPos = src.find('.');
        std::string gname = (dotPos == std::string::npos) ? src : src.substr(0, dotPos);
        if (g_useGenomeIdFilter && g_includedGenomeNames.find(gname) == g_includedGenomeNames.end()) {
            continue;
        }
        uint8_t gID = getGenomeID(gname);

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
                } else {
                    mafFile.seekg(savePos);
                }
            }
        }

        std::string cleaned = removeDashes(seq);
        args.totalLengthProcessed += cleaned.length();
        args.totalSequencesProcessed += 1;
        if (cleaned.size() < (size_t)args.k) continue;

        for (size_t pos = 0; pos + args.k <= cleaned.size(); pos++) {
            uint128_t kmerBits = 0U;
            bool valid = true;
            for (int i = 0; i < args.k; i++) {
                char c = cleaned[pos + i];
                uint128_t val;
                switch (c) {
                    case 'A': val = 0U; break;
                    case 'C': val = 1U; break;
                    case 'G': val = 2U; break;
                    case 'T': val = 3U; break;
                    default:  valid = false; break;
                }
                if (!valid) break;
                if (!cleanedQ.empty()) {
                    if (!isQValid(cleanedQ[pos + i])) {
                        valid = false;
                        break;
                    }
                }
                kmerBits = (kmerBits << 2) | val;
            }
            if (!valid) continue;

            if (g_useCanonical) {
                kmerBits = canonicalKmer128(kmerBits, args.k);
            }

            uint16_t p1 = extractP1(kmerBits, args.k);
            localBins[p1].records.push_back({kmerBits, gID});
            localCount++;

            if (localBins[p1].records.size() >= BIN_CAPACITY) {
                BinChunk chunk = partialSortAndCompact(localBins[p1], args.k);
                localBins[p1].records.clear();
                int pmIndex = (int)(chunk.p1_id % args.numPMs);
                auto &pm = (*args.pms)[pmIndex];
                {
                    std::lock_guard<std::mutex> lk(pm.pmMutex);
                    pm.pmQueue.push(std::move(chunk));
                    pm.pmCV.notify_one();
                }
            }
        }
    }
    for (auto &b : localBins) {
        if (!b.records.empty()) {
            BinChunk chunk = partialSortAndCompact(b, args.k);
            b.records.clear();
            int pmIndex = (int)(chunk.p1_id % args.numPMs);
            auto &pm = (*args.pms)[pmIndex];
            {
                std::lock_guard<std::mutex> lk(pm.pmMutex);
                pm.pmQueue.push(std::move(chunk));
                pm.pmCV.notify_one();
            }
        }
    }
    args.totalKmers->fetch_add(localCount, std::memory_order_relaxed);

    auto end = std::chrono::steady_clock::now();
    double duration = std::chrono::duration<double>(end - start).count();
    {
        std::lock_guard<std::mutex> lock(coutMutex);
        std::cout << "Reader thread " << args.readerId << " finished in " << duration << " seconds.\n";
    }

    if (args.totalSequencesProcessed > 0) {
        std::cout << "Avg sequence size processed: " 
                  << args.totalLengthProcessed / args.totalSequencesProcessed << "\n";
    }
}

// ---------------------------------------------
// Merge multiple bin_X_*.bin => bin_X_complete => final.bin
// ---------------------------------------------
struct CombinedData {
    uint8_t suffixMode;
    std::vector<uint32_t> v0;
    std::vector<uint64_t> v1;
    std::vector<uint128_t> v3;
};

struct GroupMeta {
    uint16_t p2_val;
    uint32_t sz;
};

static void processBinFiles(uint16_t p1, const std::vector<std::string> &files, const std::string &outName)
{
    bool firstFile = true;
    bool modeError = false;
    uint32_t file_k = 0;
    uint16_t file_p1_bits = 0;
    uint16_t file_p2_bits = 0;
    uint8_t suffix_mode = 0;

    std::unordered_map<uint16_t, CombinedData> combined;
    for (auto &fname : files) {
        std::ifstream ifs(fname, std::ios::binary);
        if(!ifs) continue;

        struct FileHeader {
            char     magic[4];
            uint32_t k;
            uint16_t p1_bits;
            uint16_t p2_bits;
            uint8_t  genome_id_bits;
            uint8_t  suffix_mode;
            uint8_t  reserved[2];
        } hdr;

        ifs.read((char*)&hdr, sizeof(hdr));
        if (std::strncmp(hdr.magic, "KMHD", 4) != 0) {
            continue;
        }

        uint8_t fileMode = hdr.suffix_mode;

        if (firstFile) {
            file_k        = hdr.k;
            file_p1_bits  = hdr.p1_bits;
            file_p2_bits  = hdr.p2_bits;
            suffix_mode   = fileMode;
            firstFile     = false;
        } else {
            if (hdr.k != file_k || hdr.p1_bits != file_p1_bits ||
                hdr.p2_bits != file_p2_bits || fileMode != suffix_mode)
            {
                modeError = true;
                continue;
            }
        }

        char sentinel[4];
        ifs.read(sentinel, 4);

        uint32_t groupCount = 0;
        ifs.read((char*)&groupCount, sizeof(groupCount));

        std::vector<GroupMeta> gmeta(groupCount);
        for (uint32_t i = 0; i < groupCount; i++) {
            ifs.read((char*)&gmeta[i].p2_val, sizeof(gmeta[i].p2_val));
            ifs.read((char*)&gmeta[i].sz,     sizeof(gmeta[i].sz));
        }

        ifs.read(sentinel, 4);

        for (auto &gm : gmeta) {
            if (gm.sz == 0) continue;
            auto &cd = combined[gm.p2_val];
            cd.suffixMode = suffix_mode;
            switch (suffix_mode) {
                case 0: {
                    size_t oldSize = cd.v0.size();
                    cd.v0.resize(oldSize + gm.sz);
                    ifs.read((char*)&cd.v0[oldSize], gm.sz * sizeof(uint32_t));
                } break;
                case 1: {
                    size_t oldSize = cd.v1.size();
                    cd.v1.resize(oldSize + gm.sz);
                    for (size_t idx = 0; idx < gm.sz; idx++) {
                        ifs.read((char*)&cd.v1[oldSize + idx], sizeof(uint64_t));
                    }
                } break;
                case 3: {
                    size_t oldSize = cd.v3.size();
                    cd.v3.resize(oldSize + gm.sz);
                    for (size_t idx = 0; idx < gm.sz; idx++) {
                        unsigned char bytes[16];
                        ifs.read((char*)bytes, 16);
                        uint128_t val = 0;
                        for (int b = 0; b < 16; b++) {
                            uint128_t part = (uint128_t)bytes[b];
                            val |= (part << (8ULL * b));
                        }
                        cd.v3[oldSize + idx] = val;
                    }
                } break;
            }
        }
    }

    if (modeError) {
        return;
    }

    std::ofstream ofs(outName, std::ios::binary);
    if(!ofs) return;

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
    hdr.k = file_k;
    hdr.p1_bits = file_p1_bits;
    hdr.p2_bits = file_p2_bits;
    hdr.genome_id_bits = 8;
    hdr.suffix_mode = suffix_mode;
    hdr.reserved[0] = 0; 
    hdr.reserved[1] = 0;
    ofs.write((const char*)&hdr, sizeof(hdr));
    const char *sent1 = "kmhd";
    ofs.write(sent1, 4);

    std::vector<uint16_t> p2vals;
    p2vals.reserve(combined.size());
    for (auto &kv : combined) {
        p2vals.push_back(kv.first);
    }
    std::sort(p2vals.begin(), p2vals.end());

    std::vector<uint32_t> groupSizes;
    groupSizes.reserve(p2vals.size());
    for (uint16_t v : p2vals) {
        auto &cd = combined[v];
        switch (cd.suffixMode) {
            case 0: groupSizes.push_back((uint32_t)cd.v0.size()); break;
            case 1: groupSizes.push_back((uint32_t)cd.v1.size()); break;
            case 3: groupSizes.push_back((uint32_t)cd.v3.size()); break;
        }
    }
    uint32_t groupCount = (uint32_t)p2vals.size();
    ofs.write((const char*)&groupCount, sizeof(groupCount));

    for (size_t i = 0; i < groupCount; i++) {
        uint16_t p2v = p2vals[i];
        ofs.write((const char*)&p2v, sizeof(p2v));
        uint32_t sz = groupSizes[i];
        ofs.write((const char*)&sz, sizeof(sz));
    }
    const char* sent2 = "kmhd";
    ofs.write(sent2, 4);

    // Now sort each group and write
    for (uint16_t v : p2vals) {
        auto &cd = combined[v];
        switch (cd.suffixMode) {
            case 0: {
                if (!cd.v0.empty()) {
                    LSDRadixSort32(cd.v0);
                    ofs.write((const char*)cd.v0.data(), cd.v0.size() * sizeof(uint32_t));
                }
            } break;
            case 1: {
                if (!cd.v1.empty()) {
                    LSDRadixSort64(cd.v1);
                    ofs.write((const char*)cd.v1.data(), cd.v1.size() * sizeof(uint64_t));
                }
            } break;
            case 3: {
                if (!cd.v3.empty()) {
                    LSDRadixSort128(cd.v3);
                    for (auto &val : cd.v3) {
                        unsigned char bytes[16];
                        for (int b = 0; b < 16; b++) {
                            uint64_t part = (val >> (8ULL * b)).convert_to<uint64_t>();
                            bytes[b] = static_cast<unsigned char>(part & 0xFF);
                        }
                        ofs.write((char*)bytes, 16);
                    }
                }
            } break;
        }
    }
    ofs.close();
}
static void purgeIntermediateParts(const std::string &excludeFile = "")
{
    for (auto &entry : std::filesystem::directory_iterator(g_tempFilesDir)) {
        if (!entry.is_regular_file())
            continue;
        std::string fname = entry.path().filename().string();
        // Delete files that start with "bin_" and do not contain "complete"
        // but only in g_tempFilesDir
        if (fname.rfind("bin_", 0) == 0 && fname.find("complete") == std::string::npos) {
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

static void purgeIntermediateFiles(const std::string &excludeFile = "")
{
    // Purge any file in g_tempFilesDir that starts with "bin_"
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

int main(int argc, char* argv[])
{
    // These will hold the user-provided values for thread parameters
    int userProvidedReaders = -1;
    int userProvidedPMs = -1;
    int userProvidedThreads = -1;

    std::vector<std::string> positionalArgs;

    // ------------------------------------------------------
    // PARSE COMMAND-LINE ARGUMENTS (NEW STYLE)
    //
    // Explanation of new parameters:
    //  --k <VAL>                   : The k-mer size (required)
    //  --reader_threads <VAL>      : Number of reader threads (mutually exclusive with --threads)
    //  --package_manager_threads <VAL> : Number of package manager threads (mutually exclusive with --threads)
    //  --threads <VAL>             : Total threads to be split approx as 2/3 readers, 1/3 PM (mutually exclusive with the above two)
    //  --temp_files_dir <DIR>      : Directory for temporary/intermediate files (defaults to current directory if not specified)
    //  --output_directory <DIR>    : Directory for final.bin and final.metadata (defaults to current directory if not specified)
    //  All other parameters remain as in the original usage.
    //
    // If the user omits both --reader_threads/--package_manager_threads and --threads, it will be an error.
    // ------------------------------------------------------
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
        // New named parameter for k
        else if (arg == "--k") {
            if (i + 1 < argc) {
                g_k = std::stoi(argv[++i]);
            } else {
                std::cerr << "Missing argument for --k\n";
                return 1;
            }
        }
        // New named parameter for reader threads
        else if (arg == "--reader_threads") {
            if (i + 1 < argc) {
                userProvidedReaders = std::stoi(argv[++i]);
            } else {
                std::cerr << "Missing argument for --reader_threads\n";
                return 1;
            }
        }
        // New named parameter for package manager threads
        else if (arg == "--package_manager_threads") {
            if (i + 1 < argc) {
                userProvidedPMs = std::stoi(argv[++i]);
            } else {
                std::cerr << "Missing argument for --package_manager_threads\n";
                return 1;
            }
        }
        // New named parameter for total threads (split into readers and PM)
        else if (arg == "--threads") {
            if (i + 1 < argc) {
                userProvidedThreads = std::stoi(argv[++i]);
            } else {
                std::cerr << "Missing argument for --threads\n";
                return 1;
            }
        }
        // New named parameter for temp files directory
        else if (arg == "--temp_files_dir") {
            if (i + 1 < argc) {
                g_tempFilesDir = argv[++i];
            } else {
                std::cerr << "Missing argument for --temp_files_dir\n";
                return 1;
            }
        }
        // New named parameter for output directory
        else if (arg == "--output_directory") {
            if (i + 1 < argc) {
                g_outputDirectory = argv[++i];
            } else {
                std::cerr << "Missing argument for --output_directory\n";
                return 1;
            }
        }
        else {
            // Positional args (e.g., the MAF file)
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

    // Show usage if no MAF file is given
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

    // The MAF file is the first (and only) positional argument
    std::string mafFileOriginal = positionalArgs[0];

    // Now decide how many threads for readers and PM.
    if (userProvidedThreads != -1 && (userProvidedReaders != -1 || userProvidedPMs != -1)) {
        std::cerr << "Error: --threads is mutually exclusive with --reader_threads and --package_manager_threads\n";
        return 1;
    }
    else if (userProvidedThreads != -1) {
        // Split threads: 2/3 readers, 1/3 PM
        if (userProvidedThreads <= 0) {
            std::cerr << "Error: --threads must be > 0\n";
            return 1;
        }
        // Simple approach:
        g_numReaders = (2 * userProvidedThreads) / 3; 
        if (g_numReaders < 1) g_numReaders = 1;
        g_numPMs = userProvidedThreads - g_numReaders; 
        if (g_numPMs < 1) g_numPMs = 1;
    }
    else {
        // Then user must have provided both reader_threads and package_manager_threads
        if (userProvidedReaders <= 0 || userProvidedPMs <= 0) {
            std::cerr << "Error: must provide --reader_threads <val> and --package_manager_threads <val>, or use --threads\n";
            return 1;
        }
        g_numReaders = userProvidedReaders;
        g_numPMs = userProvidedPMs;
    }

    // Check if k was set
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

    // If the user wants to purge intermediate files from a previous run, do it now
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

        // final.bin path is now in output directory with user-specified name
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

        std::vector<uint16_t> allBins;
        allBins.reserve(1 << P1_BITS);
        for (uint16_t i = 0; i < (1 << P1_BITS); i++) {
            allBins.push_back(i);
        }
        std::vector<std::string> allFiles;
        // Now we only scan g_tempFilesDir for intermediate files
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
        size_t threadCount = g_numPMs + g_numReaders;
        if (threadCount < 1) threadCount = 1;

        for (size_t i = 0; i < threadCount; i++) {
            binThreads.emplace_back([i, &outMutex, &allBins, &binMap, &nextBinIdx]() {
                auto start = std::chrono::steady_clock::now();
                while (true) {
                    uint16_t p1val = 0xFFFF;
                    {
                        std::lock_guard<std::mutex> lk(outMutex);
                        if (nextBinIdx >= allBins.size()) break;
                        p1val = allBins[nextBinIdx++];
                    }
                    auto it = binMap.find(p1val);
                    if (it == binMap.end() || it->second.empty()) {
                        continue;
                    }
                    // Construct a "complete" bin file name in temp dir
                    std::string outName = g_tempFilesDir + "/bin_" + std::to_string(p1val) + "_complete.bin";
                    processBinFiles(p1val, it->second, outName);
                    if (std::filesystem::exists(outName)) {
                        uint64_t fileSize = std::filesystem::file_size(outName);
                        g_binFileSizes[p1val] = fileSize;
                    }
                }
                auto end = std::chrono::steady_clock::now();
                double duration = std::chrono::duration<double>(end - start).count();
                {
                    std::lock_guard<std::mutex> lock(coutMutex);
                    std::cout << "Bin thread " << i << " finished in " << duration << " seconds.\n";
                }
            });
        }
        for (auto &t : binThreads) {
            t.join();
        }

        // Purge partial bin_X_*, not the "complete" ones
        std::thread asyncPurge([](){
            purgeIntermediateParts("TEST");
        });
        asyncPurge.detach();

        uint64_t totalBinSizes = 0;
        for (auto size : g_binFileSizes) {
            totalBinSizes += size;
        }
        // Calculate the header overhead.
        uint64_t overhead = 0;
        overhead += 4;               // "BINF"
        overhead += sizeof(g_k);     // g_k
        overhead += 1;               // nIDs

        for (const auto &name : g_genomeNames) {
            overhead += 2;           // length
            overhead += name.size();
        }
        overhead += 4;   // "ENDG"

        int binCount = 0;
        for (auto size : g_binFileSizes) {
            if (size > 0) {
                binCount++;
            }
        }
        overhead += binCount * (4 + 2 + 4); // Each bin chunk has overhead in final.bin

        uint64_t predictedFinalSize = totalBinSizes + overhead;
        std::cout << "Predicted final.bin size: " << predictedFinalSize << " bytes.\n";

        auto finalStart = std::chrono::steady_clock::now();

        // final.bin path is now in output directory with user-specified name
        std::string finalBinPath = g_outputDirectory + "/" + g_binaryOutputFile;
        int fd = open(finalBinPath.c_str(), O_RDWR | O_CREAT, 0666);
        if (fd < 0) {
            std::cerr << "Error creating " << finalBinPath << " via open\n";
            if (isGZ) removeFileIfExists(decompressedFile);
            return 0;
        }
        if (ftruncate(fd, predictedFinalSize) != 0) {
            std::cerr << "Error setting " << finalBinPath << " size with ftruncate\n";
            close(fd);
            if (isGZ) removeFileIfExists(decompressedFile);
            return 0;
        }

        void* map_ptr = mmap(nullptr, predictedFinalSize, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
        if (map_ptr == MAP_FAILED) {
            std::cerr << "Error: mmap failed for " << finalBinPath << "\n";
            close(fd);
            if (isGZ) removeFileIfExists(decompressedFile);
            return 0;
        }
        char* finalData = static_cast<char*>(map_ptr);
        size_t offset = 0;

        std::memcpy(finalData + offset, "BINF", 4);
        offset += 4;
        std::memcpy(finalData + offset, &g_k, sizeof(g_k));
        offset += sizeof(g_k);
        uint8_t nIDs = (uint8_t)g_genomeNames.size();
        std::memcpy(finalData + offset, &nIDs, sizeof(nIDs));
        offset += sizeof(nIDs);
        for (const auto &nm : g_genomeNames) {
            uint16_t len = (uint16_t)nm.size();
            std::memcpy(finalData + offset, &len, sizeof(len));
            offset += sizeof(len);
            std::memcpy(finalData + offset, nm.data(), nm.size());
            offset += nm.size();
        }
        std::memcpy(finalData + offset, "ENDG", 4);
        offset += 4;

        struct BinInfo {
            uint16_t p1;
            uint64_t fileSize;
            uint64_t segmentOffset;
        };

        std::vector<BinInfo> binInfos;
        binInfos.reserve(binCount);

        for (uint16_t i = 0; i < (1 << P1_BITS); i++) {
            if (g_binFileSizes[i] > 0) {
                BinInfo info;
                info.p1 = i;
                info.fileSize = g_binFileSizes[i];
                info.segmentOffset = offset;
                offset += 4 + 2 + info.fileSize + 4;
                binInfos.push_back(info);
            }
        }

        if (offset != predictedFinalSize) {
            std::cerr << "Warning: computed final size (" << predictedFinalSize 
                      << ") does not match offset (" << offset << ")\n";
        }

        std::atomic<size_t> nextIndex(0);
        auto worker = [&]() {
            while (true) {
                size_t idx = nextIndex.fetch_add(1);
                if (idx >= binInfos.size()) break;
                BinInfo info = binInfos[idx];
                size_t segOffset = info.segmentOffset;
                std::memcpy(finalData + segOffset, "BINc", 4);
                segOffset += 4;
                std::memcpy(finalData + segOffset, &info.p1, sizeof(info.p1));
                segOffset += sizeof(info.p1);
                // The "complete" file is in g_tempFilesDir
                std::string binFileName = g_tempFilesDir + "/bin_" + std::to_string(info.p1) + "_complete.bin";
                std::ifstream ifs(binFileName, std::ios::binary);
                if (ifs) {
                    size_t bytesToCopy = info.fileSize;
                    size_t bytesCopied = 0;
                    while (bytesCopied < bytesToCopy) {
                        ifs.read(finalData + segOffset + bytesCopied, bytesToCopy - bytesCopied);
                        size_t count = ifs.gcount();
                        if (count == 0) {
                            std::cerr << "Warning: unexpected EOF when reading " << binFileName << "\n";
                            break;
                        }
                        bytesCopied += count;
                    }
                    segOffset += bytesCopied;
                } else {
                    std::cerr << "Error: cannot open " << binFileName << "\n";
                }
                std::memcpy(finalData + segOffset, "BEND", 4);
            }
        };

        std::vector<std::thread> finalThreads;
        for (size_t i = 0; i < (size_t)g_numReaders; i++) {
            finalThreads.emplace_back(worker);
        }
        for (auto &t : finalThreads) {
            t.join();
        }

        if (msync(finalData, predictedFinalSize, MS_SYNC) != 0) {
            std::cerr << "Warning: msync failed\n";
        }
        munmap(finalData, predictedFinalSize);
        close(fd);

        std::cout << g_binaryOutputFile << " created using mmap at " << finalBinPath << ".\n";

        // Now also write final.metadata
        {
            std::string metaOutName = g_outputDirectory + "/final.metadata";
            std::ofstream metaOut(metaOutName);
            if (!metaOut) {
                std::cerr << "Error creating " << metaOutName << "\n";
            } else {
                for (const auto &info : binInfos) {
                    metaOut << info.p1 << " " << info.segmentOffset << "\n";
                }
            }
            metaOut.close();
        }

        auto finalEnd = std::chrono::steady_clock::now();
        double finalDuration = std::chrono::duration<double>(finalEnd - finalStart).count();
        {
            std::lock_guard<std::mutex> lock(coutMutex);
            std::cout << "Final.bin creation took " << finalDuration << " seconds.\n";
        }
        std::cout << finalBinPath << " created.\n";
    }

    auto pre_remove = std::chrono::steady_clock::now();

    if (isGZ) {
        removeFileIfExists(decompressedFile);
    }
    if (g_purgeIntermediate) {
        // We remove all "bin_" files in the temp directory
        purgeIntermediateFiles(g_binaryOutputFile);
    }

    auto post_remove = std::chrono::steady_clock::now();
    double final_remove = std::chrono::duration<double>(post_remove - pre_remove).count();
    std::cout << "Remove took" << final_remove << " seconds.\n";

    return 0;
}
