#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <thread>
#include <mutex>
#include <atomic>
#include <unordered_map>
#include <sstream>
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <chrono>
#include <google/dense_hash_map>
#include <limits>
#include <zlib.h>

// For 128-bit support
#include <boost/multiprecision/cpp_int.hpp>
using uint128_t = boost::multiprecision::uint128_t;

// For .gz decompression
#include <zlib.h>
#include <cstdio>   // for std::remove
#include <filesystem>

// -----------------------
// 5-bit encoding for AAs
// -----------------------
static const std::unordered_map<char, uint8_t> aminoAcidToCode = {
    {'A', 0}, {'R', 1}, {'N', 2}, {'D', 3}, {'C', 4},
    {'Q', 5}, {'E', 6}, {'G', 7}, {'H', 8}, {'I', 9},
    {'L', 10}, {'K', 11}, {'M', 12}, {'F', 13}, {'P', 14},
    {'S', 15}, {'T', 16}, {'W', 17}, {'Y', 18}, {'V', 19},
    {'B', 20}, {'Z', 21}, {'U', 22}, {'X', 23}, {'O', 24}
};

static int g_k = 0;
static std::string g_outputFile = "final_proteome.bin";
static int g_numThreads = 1;

// Detect if file ends with ".gz"
bool endsWithGZ(const std::string &filename) {
    return (filename.size() >= 3 &&
            filename.compare(filename.size() - 3, 3, ".gz") == 0);
}

std::string decompressGZ(const std::string &gzFilename)
{
    // Choose a random/unique name
    auto now = std::chrono::steady_clock::now().time_since_epoch().count();
    std::string tempName = "temp_unzipped_" + std::to_string(now) + ".maf";

    gzFile in = gzopen(gzFilename.c_str(), "rb");
    if (!in) {
        std::cerr << "Error: cannot open gzip file " << gzFilename << "\n";
        return {};
    }
    FILE *out = std::fopen(tempName.c_str(), "wb");
    if (!out) {
        std::cerr << "Error: cannot create temp file " << tempName << "\n";
        gzclose(in);
        return {};
    }
    const size_t BUF_SIZE = 65536;
    std::vector<char> buffer(BUF_SIZE);
    int bytesRead = 0;
    while ((bytesRead = gzread(in, buffer.data(), (unsigned)buffer.size())) > 0) {
        std::fwrite(buffer.data(), 1, bytesRead, out);
    }
    gzclose(in);
    std::fclose(out);

    return tempName;
}


// Remove dashes
static std::string removeDashes(const std::string &alignedSeq) {
    std::string result;
    result.reserve(alignedSeq.size());
    for (char c : alignedSeq) {
        if (c != '-') result.push_back(c);
    }
    return result;
}

// Genome ID (0..255) assignment
static std::mutex g_genomeMutex;
static std::unordered_map<std::string, uint8_t> g_genomeToId;
static std::atomic<uint8_t> g_nextId(0);

static uint8_t getGenomeID(const std::string &name) {
    std::lock_guard<std::mutex> lock(g_genomeMutex);
    auto it = g_genomeToId.find(name);
    if (it != g_genomeToId.end()) {
        return it->second;
    }
    if (g_nextId.load() == 255) {
        std::cerr << "Error: more than 256 distinct proteome IDs.\n";
        std::exit(1);
    }
    uint8_t newId = g_nextId.fetch_add(1);
    g_genomeToId[name] = newId;
    return newId;
}

// 5-bit encoding for KeyT
template <typename KeyT>
KeyT encodeAAKmer(const std::string &seq, size_t pos, int k) {
    KeyT val = 0;
    for (int i = 0; i < k; i++) {
        char aa = std::toupper(seq[pos + i]);
        auto it = aminoAcidToCode.find(aa);
        if (it == aminoAcidToCode.end()) {
            return (std::numeric_limits<KeyT>::max)(); // sentinel => invalid
        }
        val = (val << 5) | it->second;  // shift left 5 bits
    }
    return val;
}

// We'll store counts in an array of 256
using GenomeCountArray = std::array<uint32_t, 256>;

// Thread data
template <typename KeyT>
struct ThreadData {
    google::dense_hash_map<KeyT, GenomeCountArray> localMap;
    uint64_t localCount = 0;
    ThreadData() {
        // set empty/deleted
        if constexpr (std::is_same<KeyT, uint32_t>::value) {
            localMap.set_empty_key(std::numeric_limits<uint32_t>::max());
            localMap.set_deleted_key(std::numeric_limits<uint32_t>::max()-1);
        }
        else if constexpr (std::is_same<KeyT, uint64_t>::value) {
            localMap.set_empty_key(std::numeric_limits<uint64_t>::max());
            localMap.set_deleted_key(std::numeric_limits<uint64_t>::max()-1);
        }
        else {
            uint128_t emptyVal  = (uint128_t(1)<<127) - 1;
            uint128_t deleteVal = (uint128_t(1)<<127) - 2;
            localMap.set_empty_key(emptyVal);
            localMap.set_deleted_key(deleteVal);
        }
    }
};

template <typename KeyT>
struct ReaderArgs {
    std::string mafFile;
    std::streampos startPos;
    std::streampos endPos;
    int k;
    ThreadData<KeyT> *tdata;
    int tid;
};

static void findChunkBoundaries(const std::string &filename,
                                int nThreads,
                                std::vector<std::streampos> &chunks)
{
    std::ifstream in(filename);
    if(!in) {
        std::cerr << "Error: cannot open " << filename << "\n";
        std::exit(1);
    }
    in.seekg(0, std::ios::end);
    std::streampos fsize = in.tellg();
    in.seekg(0, std::ios::beg);

    chunks.resize(nThreads+1);
    for (int i=0; i<=nThreads; i++) {
        chunks[i] = fsize * i / nThreads;
    }
}

// The thread function
template <typename KeyT>
void threadFunc(ReaderArgs<KeyT> args)
{
    std::ifstream maf(args.mafFile);
    maf.seekg(args.startPos);
    std::string line;
    while ((args.endPos < 0 || maf.tellg() < args.endPos) && std::getline(maf, line)) {
        if (line.empty() || line[0] != 's') {
            continue;
        }
        // parse
        // s src start size strand srcSize sequence
        std::istringstream iss(line);
        std::string sTag, src, seq;
        int ignoreInt;
        char ignoreChar;
        iss >> sTag >> src >> ignoreInt >> ignoreInt >> ignoreChar >> ignoreInt >> seq;
        size_t dotPos = src.find('.');
        std::string genome = (dotPos == std::string::npos) ? src : src.substr(0, dotPos);
        uint8_t gid = getGenomeID(genome);
        // remove dashes
        std::string cleaned = removeDashes(seq);
        if (cleaned.size() < (size_t)args.k) continue;

        auto &mapRef = args.tdata->localMap;
        uint64_t &countRef = args.tdata->localCount;

        for (size_t i = 0; i + args.k <= cleaned.size(); i++) {
            KeyT enc = encodeAAKmer<KeyT>(cleaned, i, args.k);
            if (enc == (std::numeric_limits<KeyT>::max)()) {
                continue;
            }
            auto it = mapRef.find(enc);
            if (it == mapRef.end()) {
                GenomeCountArray arr{};
                arr[gid]++;
                mapRef.insert({enc, arr});
            } else {
                it->second[gid]++;
            }
            countRef++;
        }
    }
}

// Merge
template <typename KeyT>
void mergeLocal(
    google::dense_hash_map<KeyT, GenomeCountArray> &globalMap,
    const google::dense_hash_map<KeyT, GenomeCountArray> &localMap)
{
    for (auto &kv : localMap) {
        auto it = globalMap.find(kv.first);
        if (it == globalMap.end()) {
            globalMap.insert({kv.first, kv.second});
        } else {
            for (int i=0; i<256; i++) {
                it->second[i] += kv.second[i];
            }
        }
    }
}

// Write final output
// 4 bytes: "PRTM"
// 4 bytes: k
// 8 bytes: number_of_kmers
// For each k-mer => key(4,8,16) + array[256] (4 bytes each)
template <typename KeyT>
void writeFinal(const std::string &filename,
                const google::dense_hash_map<KeyT, GenomeCountArray> &globalMap)
{
    std::ofstream out(filename, std::ios::binary);
    if (!out) {
        std::cerr << "Error writing " << filename << "\n";
        return;
    }
    out.write("PRTM", 4);
    uint32_t k32 = (uint32_t)g_k;
    out.write((char*)&k32, sizeof(k32));
    uint64_t nKmers = globalMap.size();
    out.write((char*)&nKmers, sizeof(nKmers));
    for (auto &kv : globalMap) {
        if constexpr (std::is_same<KeyT, uint32_t>::value) {
            uint32_t val = kv.first;
            out.write((char*)&val, 4);
        }
        else if constexpr (std::is_same<KeyT, uint64_t>::value) {
            uint64_t val = kv.first;
            out.write((char*)&val, 8);
        }
        else {
            // 128 bits
            uint128_t val = kv.first;
            unsigned char bytes[16];
            for (int i = 0; i < 16; i++) {
                uint64_t part = (val >> (8ULL * i)).convert_to<uint64_t>();
                bytes[i] = static_cast<unsigned char>(part & 0xFF);
            }
            out.write((char*)bytes, 16);
        }
        // write 256 x uint32_t
        out.write((char*)kv.second.data(), 256 * sizeof(uint32_t));
    }
    out.close();
    std::cout << "Wrote " << nKmers << " k-mers to " << filename << "\n";
}

// Main
int main(int argc, char *argv[])
{
    if (argc < 3) {
        std::cerr << "Usage:\n"
                  << "  " << argv[0] << " --k <kmer_size> [--threads N]\n"
                  << "       [--binary_file_output <file>] <MAFFile>\n";
        return 1;
    }

    std::string mafFile;
    for (int i=1; i<argc; i++) {
        std::string arg(argv[i]);
        if (arg == "--k") {
            if (++i < argc) {
                g_k = std::stoi(argv[i]);
            } else {
                std::cerr << "Error: --k requires an integer\n";
                return 1;
            }
        }
        else if (arg == "--threads") {
            if (++i < argc) {
                g_numThreads = std::stoi(argv[i]);
                if (g_numThreads < 1) g_numThreads = 1;
            } else {
                std::cerr << "Error: --threads requires an integer\n";
                return 1;
            }
        }
        else if (arg == "--binary_file_output") {
            if (++i < argc) {
                g_outputFile = argv[i];
            } else {
                std::cerr << "Error: --binary_file_output requires a filename\n";
                return 1;
            }
        }
        else {
            // presumably MAF
            mafFile = arg;
        }
    }

    if (g_k < 1 || g_k > 25) {
        std::cerr << "Error: for proteome k‑mers, k must be 1..25.\n";
        return 1;
    }
    if (mafFile.empty()) {
        std::cerr << "Error: no MAF file specified.\n";
        return 1;
    }

    // Check if gz
    std::string actualFile = mafFile;
    bool gz = endsWithGZ(mafFile);
    std::string tmpFile;
    if (gz) {
        tmpFile = decompressGZ(mafFile);
        if (tmpFile.empty()) {
            std::cerr << "Error: decompression failed.\n";
            return 1;
        }
        actualFile = tmpFile;
        std::cout << "Decompressed " << mafFile << " -> " << tmpFile << "\n";
    }

    // chunk
    std::vector<std::streampos> chunkBound;
    findChunkBoundaries(actualFile, g_numThreads, chunkBound);

    // pick KeyT
    // k <= 6 => 32 bits
    // k <=12 => 64 bits
    // k <=25 => 128 bits
    if (g_k <= 6) {
        using KeyType = uint32_t;
        std::vector<ThreadData<KeyType>> tdata(g_numThreads);
        std::vector<std::thread> threads;
        for (int i=0; i<g_numThreads; i++) {
            ReaderArgs<KeyType> rargs{
                actualFile, chunkBound[i], chunkBound[i+1], g_k,
                &tdata[i], i
            };
            threads.emplace_back(threadFunc<KeyType>, rargs);
        }
        for (auto &t : threads) t.join();

        // Merge
        google::dense_hash_map<KeyType, GenomeCountArray> globalMap;
        globalMap.set_empty_key(std::numeric_limits<KeyType>::max());
        globalMap.set_deleted_key(std::numeric_limits<KeyType>::max()-1);

        uint64_t totalRaw = 0;
        for (auto &td : tdata) {
            mergeLocal<KeyType>(globalMap, td.localMap);
            totalRaw += td.localCount;
        }
        std::cout << "Total raw k‑mers processed: " << totalRaw
                  << "; unique: " << globalMap.size() << "\n";

        writeFinal<KeyType>(g_outputFile, globalMap);
    }
    else if (g_k <= 12) {
        using KeyType = uint64_t;
        std::vector<ThreadData<KeyType>> tdata(g_numThreads);
        std::vector<std::thread> threads;
        for (int i=0; i<g_numThreads; i++) {
            ReaderArgs<KeyType> rargs{
                actualFile, chunkBound[i], chunkBound[i+1], g_k,
                &tdata[i], i
            };
            threads.emplace_back(threadFunc<KeyType>, rargs);
        }
        for (auto &t : threads) t.join();

        google::dense_hash_map<KeyType, GenomeCountArray> globalMap;
        globalMap.set_empty_key(std::numeric_limits<KeyType>::max());
        globalMap.set_deleted_key(std::numeric_limits<KeyType>::max()-1);

        uint64_t totalRaw = 0;
        for (auto &td : tdata) {
            mergeLocal<KeyType>(globalMap, td.localMap);
            totalRaw += td.localCount;
        }
        std::cout << "Total raw k‑mers processed: " << totalRaw
                  << "; unique: " << globalMap.size() << "\n";

        writeFinal<KeyType>(g_outputFile, globalMap);
    }
    else {
        // up to k=25 => 128 bits
        using KeyType = uint128_t;
        std::vector<ThreadData<KeyType>> tdata(g_numThreads);
        std::vector<std::thread> threads;
        for (int i=0; i<g_numThreads; i++) {
            ReaderArgs<KeyType> rargs{
                actualFile, chunkBound[i], chunkBound[i+1], g_k,
                &tdata[i], i
            };
            threads.emplace_back(threadFunc<KeyType>, rargs);
        }
        for (auto &t : threads) t.join();

        google::dense_hash_map<KeyType, GenomeCountArray> globalMap;
        {
            uint128_t emptyVal  = (uint128_t(1)<<127) - 1;
            uint128_t deleteVal = (uint128_t(1)<<127) - 2;
            globalMap.set_empty_key(emptyVal);
            globalMap.set_deleted_key(deleteVal);
        }

        uint64_t totalRaw = 0;
        for (auto &td : tdata) {
            mergeLocal<KeyType>(globalMap, td.localMap);
            totalRaw += td.localCount;
        }
        std::cout << "Total raw k‑mers processed: " << totalRaw
                  << "; unique: " << globalMap.size() << "\n";

        writeFinal<KeyType>(g_outputFile, globalMap);
    }

    // remove temp if needed
    if (gz && !tmpFile.empty()) {
        std::error_code ec;
        std::filesystem::remove(tmpFile, ec);
    }
    {
    std::ofstream metaOut("proteomes.metadata");
    if (!metaOut) {
        std::cerr << "Error: cannot open proteomes.metadata for writing.\n";
        return 1;
    }
    // Write header if desired (e.g. "Proteome\tID")
    metaOut << "Proteome\tID\n";
    for (const auto &entry : g_genomeToId) {
        metaOut << entry.first << "\t" << static_cast<int>(entry.second) << "\n";
    }
    metaOut.close();
    std::cout << "proteomes.metadata created.\n";
}


    return 0;
} 