// kmer_counter.cpp

#include <sparsehash/dense_hash_map>
#include <sparsehash/dense_hash_set>
#include <unordered_map>
#include <vector>
#include <string>
#include <cstdint>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>     // For exit()
#include <sys/stat.h>  // For mkdir
#include <sys/types.h> // For mkdir
#include <errno.h>     // For errno
#include <cstring>     // For strerror
#include <chrono>      // For timing
#include <thread>      // For multithreading
#include <mutex>       // For mutex
#include <algorithm>   // For min
#include "concurrentqueue.h" // For moodycamel::ConcurrentQueue

using google::dense_hash_map;
using google::dense_hash_set;

// Custom hash functions for uint32_t and uint64_t keys
struct UInt32Hasher {
    std::size_t operator()(uint32_t key) const noexcept {
        return std::hash<uint32_t>{}(key);
    }
};

struct UInt64Hasher {
    std::size_t operator()(uint64_t key) const noexcept {
        // Custom hash function for better distribution
        key = (~key) + (key << 21);
        key = key ^ (key >> 24);
        key = (key + (key << 3)) + (key << 8);
        key = key ^ (key >> 14);
        key = (key + (key << 2)) + (key << 4);
        key = key ^ (key >> 28);
        key = key + (key << 31);
        return static_cast<std::size_t>(key);
    }
};

// Hash function selector based on KmerCodeType
template<typename T>
struct KmerCodeHasher;
template<>
struct KmerCodeHasher<uint32_t> : public UInt32Hasher {};
template<>
struct KmerCodeHasher<uint64_t> : public UInt64Hasher {};

// Custom struct to store genome counts in KmerGroups and KmerCounts
struct GenomeInfo {
    dense_hash_set<uint8_t> genomesWithCountOne;
    dense_hash_map<uint8_t, uint16_t> genomesWithMultipleCounts;

    GenomeInfo() {
        genomesWithCountOne.set_empty_key(255);
        genomesWithCountOne.set_deleted_key(254);
        genomesWithMultipleCounts.set_empty_key(253);
        genomesWithMultipleCounts.set_deleted_key(252);
    }

    void increment(uint8_t genomeId) {
        if (genomesWithCountOne.find(genomeId) != genomesWithCountOne.end()) {
            // genomeId was already in genomesWithCountOne
            genomesWithCountOne.erase(genomeId);
            genomesWithMultipleCounts[genomeId] = 2;
        } else if (genomesWithMultipleCounts.find(genomeId) != genomesWithMultipleCounts.end()) {
            // genomeId is in genomesWithMultipleCounts
            genomesWithMultipleCounts[genomeId]++;
        } else {
            // genomeId is new
            genomesWithCountOne.insert(genomeId);
        }
    }

    // Function to get all genome IDs and counts
    std::vector<std::pair<uint8_t, uint16_t>> getGenomeIdCounts() const {
        std::vector<std::pair<uint8_t, uint16_t>> genomeIdCounts;
        genomeIdCounts.reserve(genomesWithCountOne.size() + genomesWithMultipleCounts.size());

        for (uint8_t genomeId : genomesWithCountOne) {
            genomeIdCounts.emplace_back(genomeId, 1);
        }

        for (const auto& entry : genomesWithMultipleCounts) {
            genomeIdCounts.emplace_back(entry.first, entry.second);
        }

        return genomeIdCounts;
    }
};

template<typename KmerCodeType>
class KmerCounter {
public:
    // ID mappings
    std::unordered_map<std::string, uint8_t> idToIndex;
    std::vector<std::string> indexToId;
    std::mutex idMappingMutex; // Mutex for thread-safe genome ID mapping

    // K-mer counts
    typedef dense_hash_map<KmerCodeType, GenomeInfo, KmerCodeHasher<KmerCodeType>> KmerCountsMap;
    std::vector<KmerCountsMap> kmerCountsPartitions;

    // Timing variables
    std::chrono::duration<double> totalProcessSequenceTime = std::chrono::duration<double>::zero();
    std::chrono::duration<double> totalMergeTime = std::chrono::duration<double>::zero();
    std::chrono::duration<double> totalWritingTime = std::chrono::duration<double>::zero();

    // Counters
    size_t totalKMersProcessed = 0;
    size_t totalSequencesProcessed = 0;
    size_t totalAlignmentBlocks = 0;

    KmerCounter(int k, int numConsumers, bool useCanonicalKmers)
        : k(k), numConsumers(numConsumers), useCanonicalKmers(useCanonicalKmers) {
        // Initialize kmerCountsPartitions
        kmerCountsPartitions.resize(numConsumers);
        for (auto& partition : kmerCountsPartitions) {
            partition.set_empty_key(static_cast<KmerCodeType>(-1));
            partition.set_deleted_key(static_cast<KmerCodeType>(-2));
        }

        producersDone = false;
    }

    // Producer thread function to process its assigned chunk
    void processChunk(const std::string& filename, std::streampos startPos, std::streampos endPos,
                  std::vector<moodycamel::ConcurrentQueue<KmerCountsMap>>& kmerGroupQueues,
                  int threadId) {
        std::cout << "[Producer " << threadId << "] Started processing its chunk." << std::endl;

        auto chunkStartTime = std::chrono::high_resolution_clock::now();
        size_t localTotalKMersProcessed = 0;
        size_t localTotalSequencesProcessed = 0;
        size_t localTotalAlignmentBlocks = 0;

        std::ifstream mafFile(filename);
        if (!mafFile) {
            std::cerr << "Error opening MAF file: " << filename << std::endl;
            exit(1);
        }

        mafFile.seekg(startPos);

        std::string line;
        bool inAlignmentBlock = false;
        std::vector<std::pair<std::string, std::string>> sequences;

        while ((mafFile.tellg() < endPos || endPos == std::streampos(-1)) && std::getline(mafFile, line)) {
            if (line.empty() || line[0] == '#') {
                continue;
            }

            if (line[0] == 'a') {
                if (!sequences.empty()) {
                    // Process the alignment block
                    localTotalAlignmentBlocks++;

                    // Initialize per-partition KmerGroups
                    std::vector<KmerCountsMap> kmerGroups(numConsumers);
                    for (auto& kg : kmerGroups) {
                        kg.set_empty_key(static_cast<KmerCodeType>(-1));
                        kg.set_deleted_key(static_cast<KmerCodeType>(-2));
                    }

                    size_t numSequences = sequences.size();
                    localTotalSequencesProcessed += numSequences;

                    // Process each sequence individually
                    for (const auto& seqPair : sequences) {
                        const std::string& id = seqPair.first;
                        const std::string& alignedSeq = seqPair.second;

                        uint8_t idIndex;

                        // Attempt to find the genome ID without locking (not thread-safe)
                        auto it = idToIndex.find(id);
                        if (it != idToIndex.end()) {
                            // Genome ID found, proceed without locking
                            idIndex = it->second;
                        } else {
                            // Lock and double-check
                            {
                                std::unique_lock<std::mutex> lock(idMappingMutex);
                                // Double-check if another thread added the ID
                                it = idToIndex.find(id);
                                if (it == idToIndex.end()) {
                                    if (idToIndex.size() >= 254) {
                                        std::cerr << "Error: Exceeded maximum number of genome IDs (254)." << std::endl;
                                        exit(1);
                                    }
                                    idIndex = idToIndex.size();
                                    idToIndex[id] = idIndex;
                                    indexToId.push_back(id);
                                } else {
                                    idIndex = it->second;
                                }
                                // Lock is released when out of scope
                            }
                        }

                        size_t seqLength = alignedSeq.size();
                        KmerCodeType kmerCode = 0;
                        size_t validBases = 0;
                        KmerCodeType mask = (k >= sizeof(KmerCodeType) * 4) ? static_cast<KmerCodeType>(-1)
                                                                            : ((static_cast<KmerCodeType>(1) << (2 * k)) - 1);

                        for (size_t i = 0; i < seqLength; ++i) {
                            char nucleotide = alignedSeq[i];
                            KmerCodeType code = 0;

                            switch (nucleotide) {
                                case 'A': case 'a': code = 0b00; break;
                                case 'C': case 'c': code = 0b01; break;
                                case 'G': case 'g': code = 0b10; break;
                                case 'T': case 't': code = 0b11; break;
                                case '-':
                                    continue;
                                default:
                                    validBases = 0;
                                    kmerCode = 0;
                                    continue;
                            }

                            kmerCode = ((kmerCode << 2) | code) & mask;
                            validBases++;

                            if (validBases >= static_cast<size_t>(k)) {
                                KmerCodeType finalKmerCode = kmerCode;
                                if (useCanonicalKmers) {
                                    KmerCodeType revCompCode = reverseComplement(kmerCode, k);
                                    finalKmerCode = std::min(kmerCode, revCompCode);
                                }
                                int partitionIndex = finalKmerCode % numConsumers;
                                auto& genomeInfo = kmerGroups[partitionIndex][finalKmerCode];
                                genomeInfo.increment(idIndex);
                                localTotalKMersProcessed++;
                            }
                        }
                    }

                    // Push per-partition kmerGroups into the respective queues
                    for (int i = 0; i < numConsumers; ++i) {
                        if (!kmerGroups[i].empty()) {
                            kmerGroupQueues[i].enqueue(std::move(kmerGroups[i]));
                        }
                    }

                    sequences.clear();

                    // Logging progress every 1000 blocks
                    if (localTotalAlignmentBlocks % 1000 == 0) {
                        std::cout << "[Producer " << threadId << "] Processed " << localTotalAlignmentBlocks << " alignment blocks." << std::endl;
                    }
                }
                inAlignmentBlock = true;
                continue;
            }

            if (inAlignmentBlock && line[0] == 's') {
                std::istringstream iss(line);
                std::string tag;
                std::string src;
                int start, size;
                char strand;
                int srcSize;
                std::string text;

                iss >> tag >> src >> start >> size >> strand >> srcSize >> text;

                // Extract only the part of src before the first period
                src = src.substr(0, src.find('.'));
                text = text.substr(0, text.find_first_of(" \t\n"));

                sequences.emplace_back(src, text);
            } else if (line[0] != 's') {
                inAlignmentBlock = false;
            }
        }

        // Process any remaining sequences
        if (!sequences.empty()) {
            // Process the alignment block
            localTotalAlignmentBlocks++;

            // Initialize per-partition KmerGroups
            std::vector<KmerCountsMap> kmerGroups(numConsumers);
            for (auto& kg : kmerGroups) {
                kg.set_empty_key(static_cast<KmerCodeType>(-1));
                kg.set_deleted_key(static_cast<KmerCodeType>(-2));
            }

            size_t numSequences = sequences.size();
            localTotalSequencesProcessed += numSequences;

            // Process each sequence individually
            for (const auto& seqPair : sequences) {
                const std::string& id = seqPair.first;
                const std::string& alignedSeq = seqPair.second;

                uint8_t idIndex;

                // Attempt to find the genome ID without locking (not thread-safe)
                auto it = idToIndex.find(id);
                if (it != idToIndex.end()) {
                    // Genome ID found, proceed without locking
                    idIndex = it->second;
                } else {
                    // Lock and double-check
                    {
                        std::unique_lock<std::mutex> lock(idMappingMutex);
                        // Double-check if another thread added the ID
                        it = idToIndex.find(id);
                        if (it == idToIndex.end()) {
                            if (idToIndex.size() >= 254) {
                                std::cerr << "Error: Exceeded maximum number of genome IDs (254)." << std::endl;
                                exit(1);
                            }
                            idIndex = idToIndex.size();
                            idToIndex[id] = idIndex;
                            indexToId.push_back(id);
                        } else {
                            idIndex = it->second;
                        }
                        // Lock is released when out of scope
                    }
                }

                size_t seqLength = alignedSeq.size();
                KmerCodeType kmerCode = 0;
                size_t validBases = 0;
                KmerCodeType mask = (k >= sizeof(KmerCodeType) * 4) ? static_cast<KmerCodeType>(-1)
                                                                    : ((static_cast<KmerCodeType>(1) << (2 * k)) - 1);

                for (size_t i = 0; i < seqLength; ++i) {
                    char nucleotide = alignedSeq[i];
                    KmerCodeType code = 0;

                    switch (nucleotide) {
                        case 'A': case 'a': code = 0b00; break;
                        case 'C': case 'c': code = 0b01; break;
                        case 'G': case 'g': code = 0b10; break;
                        case 'T': case 't': code = 0b11; break;
                        case '-':
                            continue;
                        default:
                            validBases = 0;
                            kmerCode = 0;
                            continue;
                    }

                    kmerCode = ((kmerCode << 2) | code) & mask;
                    validBases++;

                    if (validBases >= static_cast<size_t>(k)) {
                        KmerCodeType finalKmerCode = kmerCode;
                        if (useCanonicalKmers) {
                            KmerCodeType revCompCode = reverseComplement(kmerCode, k);
                            finalKmerCode = std::min(kmerCode, revCompCode);
                        }
                        int partitionIndex = finalKmerCode % numConsumers;
                        auto& genomeInfo = kmerGroups[partitionIndex][finalKmerCode];
                        genomeInfo.increment(idIndex);
                        localTotalKMersProcessed++;
                    }
                }
            }

            // Push per-partition kmerGroups into the respective queues
            for (int i = 0; i < numConsumers; ++i) {
                if (!kmerGroups[i].empty()) {
                    kmerGroupQueues[i].enqueue(std::move(kmerGroups[i]));
                }
            }

            sequences.clear();

            // Logging progress
            if (localTotalAlignmentBlocks % 1000 == 0) {
                std::cout << "[Producer " << threadId << "] Processed " << localTotalAlignmentBlocks << " alignment blocks." << std::endl;
            }
        }

        mafFile.close();

        // Update global counters
        totalKMersProcessed += localTotalKMersProcessed;
        totalSequencesProcessed += localTotalSequencesProcessed;
        totalAlignmentBlocks += localTotalAlignmentBlocks;

        auto chunkEndTime = std::chrono::high_resolution_clock::now();
        auto chunkDuration = std::chrono::duration_cast<std::chrono::seconds>(chunkEndTime - chunkStartTime).count();

        std::cout << "[Producer " << threadId << "] Finished processing its chunk. Time taken: " << chunkDuration << " seconds." << std::endl;
        std::cout << "[Producer " << threadId << "] Total kmers processed: " << localTotalKMersProcessed << std::endl;
    }


    // Consumer thread function remains unchanged
    void consumeKmerGroups(moodycamel::ConcurrentQueue<KmerCountsMap>& myQueue, int threadId) {
        std::cout << "[Consumer " << threadId << "] Started." << std::endl;
        KmerCountsMap kmerGroups;
        size_t localMergeCount = 0;

        while (true) {
            while (myQueue.try_dequeue(kmerGroups)) {
                // Merge kmerGroups into kmerCountsPartitions[threadId]
                auto mergeStart = std::chrono::high_resolution_clock::now();

                for (const auto& entry : kmerGroups) {
                    KmerCodeType kmerCode = entry.first;
                    const GenomeInfo& blockGenomeInfo = entry.second;

                    auto& genomeInfo = kmerCountsPartitions[threadId][kmerCode];

                    // Merge genomesWithCountOne
                    for (uint8_t genomeId : blockGenomeInfo.genomesWithCountOne) {
                        if (genomeInfo.genomesWithCountOne.find(genomeId) != genomeInfo.genomesWithCountOne.end()) {
                            // Move to genomesWithMultipleCounts
                            genomeInfo.genomesWithCountOne.erase(genomeId);
                            genomeInfo.genomesWithMultipleCounts[genomeId] = 2;
                        } else if (genomeInfo.genomesWithMultipleCounts.find(genomeId) != genomeInfo.genomesWithMultipleCounts.end()) {
                            genomeInfo.genomesWithMultipleCounts[genomeId]++;
                        } else {
                            genomeInfo.genomesWithCountOne.insert(genomeId);
                        }
                    }

                    // Merge genomesWithMultipleCounts
                    for (const auto& countEntry : blockGenomeInfo.genomesWithMultipleCounts) {
                        uint8_t genomeId = countEntry.first;
                        uint16_t count = countEntry.second;

                        if (genomeInfo.genomesWithCountOne.find(genomeId) != genomeInfo.genomesWithCountOne.end()) {
                            genomeInfo.genomesWithCountOne.erase(genomeId);
                            genomeInfo.genomesWithMultipleCounts[genomeId] = count + 1;
                        } else {
                            genomeInfo.genomesWithMultipleCounts[genomeId] += count;
                        }
                    }
                }

                auto mergeEnd = std::chrono::high_resolution_clock::now();
                totalMergeTime += mergeEnd - mergeStart;
                localMergeCount++;

                // Logging progress every 1000 merges
                if (localMergeCount % 1000 == 0) {
                    std::cout << "[Consumer " << threadId << "] Merged " << localMergeCount << " KmerGroups." << std::endl;
                }
            }

            // Exit the loop if producers are done and the queue is empty
            if (producersDone && myQueue.size_approx() == 0) {
                break;
            }

            // Sleep or yield to prevent busy-waiting
            std::this_thread::yield();
        }
        std::cout << "[Consumer " << threadId << "] Finished merging. Total KmerGroups merged: " << localMergeCount << std::endl;
    }

    // Function to signal that producer threads are done
    void setProducersDone() {
        producersDone = true;
    }

    // Function to write output files
    void writeOutputFiles(int numThreads, bool singleFileOutput) {
        auto start = std::chrono::high_resolution_clock::now();

        std::cout << "Writing output files using " << numThreads << " threads..." << std::endl;

        // Create the results_counter directory
        const char* dirName = "results_counter";
        struct stat st = {0};
        if (stat(dirName, &st) == -1) {
            if (mkdir(dirName, 0755) != 0) {
                std::cerr << "Error creating directory " << dirName << ": " << strerror(errno) << std::endl;
                exit(1);
            }
        }

        size_t totalKmersWritten = 0;
        size_t uniqueKmers = 0;

        // Mutex for updating shared counters
        std::mutex countMutex;

        if (singleFileOutput) {
            // Write to a single file
            std::string filename = std::string(dirName) + "/kmer_counts.txt";
            std::ofstream outFile(filename);
            if (!outFile) {
                std::cerr << "Error opening output file: " << filename << std::endl;
                exit(1);
            }

            // Mutex for thread-safe writing to the single file
            std::mutex writeMutex;

            // Create threads for writing
            std::vector<std::thread> writerThreads;

            // Distribute k-mers among threads
            size_t totalPartitions = kmerCountsPartitions.size();
            size_t partitionsPerThread = (totalPartitions + numThreads - 1) / numThreads;

            for (int threadId = 0; threadId < numThreads; ++threadId) {
                size_t startPartition = threadId * partitionsPerThread;
                size_t endPartition = std::min(startPartition + partitionsPerThread, totalPartitions);

                writerThreads.emplace_back([this, startPartition, endPartition, &outFile, &writeMutex, &totalKmersWritten, &uniqueKmers, &countMutex]() {
                    size_t localKmersWritten = 0;
                    size_t localUniqueKmers = 0;

                    for (size_t partitionId = startPartition; partitionId < endPartition; ++partitionId) {
                        auto& kmerCounts = kmerCountsPartitions[partitionId];

                        for (const auto& kmerEntry : kmerCounts) {
                            KmerCodeType kmerCode = kmerEntry.first;
                            const GenomeInfo& genomeInfo = kmerEntry.second;

                            // Decode k-mer for output
                            std::string kmerStr = decodeKmer(kmerCode);

                            // Collect genome IDs and counts
                            std::vector<std::pair<uint8_t, uint16_t>> genomeIdCounts = genomeInfo.getGenomeIdCounts();

                            // Convert genome IDs to names and prepare output line
                            std::ostringstream oss;
                            oss << kmerStr;

                            for (const auto& idCountPair : genomeIdCounts) {
                                uint8_t idIndex = idCountPair.first;
                                uint16_t count = idCountPair.second;
                                std::string genomeName = indexToId[idIndex];
                                oss << " " << genomeName << ": " << count;
                            }

                            // Write to output file
                            {
                                std::lock_guard<std::mutex> lock(writeMutex);
                                outFile << oss.str() << "\n";
                            }

                            localKmersWritten++;
                            localUniqueKmers++;
                        }
                    }

                    // Update total kmers written
                    {
                        std::lock_guard<std::mutex> lock(countMutex);
                        totalKmersWritten += localKmersWritten;
                        uniqueKmers += localUniqueKmers;
                    }
                });
            }

            // Wait for writer threads to finish
            for (auto& t : writerThreads) {
                t.join();
            }

            outFile.close();
            std::cout << "Total k-mers written to output file: " << totalKmersWritten << "\n";
        } else {
            // Distribute genomes among threads
            size_t numGenomes = indexToId.size();
            size_t genomesPerThread = (numGenomes + numThreads - 1) / numThreads;

            std::vector<std::thread> writerThreads;

            for (int threadId = 0; threadId < numThreads; ++threadId) {
                size_t startGenome = threadId * genomesPerThread;
                size_t endGenome = std::min(startGenome + genomesPerThread, numGenomes);

                writerThreads.emplace_back([this, startGenome, endGenome, dirName, threadId, &totalKmersWritten, &uniqueKmers, &countMutex]() {
                    std::cout << "[Writer " << threadId << "] Started writing files for genomes " << startGenome << " to " << endGenome - 1 << std::endl;

                    // Open output files for assigned genomes
                    std::unordered_map<uint8_t, std::ofstream> outputFiles;
                    for (size_t i = startGenome; i < endGenome; ++i) {
                        std::string filename = std::string(dirName) + "/" + indexToId[i] + "_kmer_counts.txt";
                        outputFiles[i].open(filename);
                        if (!outputFiles[i]) {
                            std::cerr << "[Writer " << threadId << "] Error opening output file: " << filename << std::endl;
                            exit(1);
                        }
                    }

                    size_t localKmersWritten = 0;
                    size_t localUniqueKmers = 0;

                    // Iterate over partitions and write kmers for assigned genomes
                    for (int partitionId = 0; partitionId < numConsumers; ++partitionId) {
                        auto& kmerCounts = kmerCountsPartitions[partitionId];

                        for (const auto& kmerEntry : kmerCounts) {
                            KmerCodeType kmerCode = kmerEntry.first;
                            const GenomeInfo& genomeInfo = kmerEntry.second;

                            // Decode k-mer for output
                            std::string kmerStr = decodeKmer(kmerCode);

                            // Write genomes with count one
                            for (uint8_t idIndex : genomeInfo.genomesWithCountOne) {
                                if (idIndex >= startGenome && idIndex < endGenome) {
                                    // Write to the corresponding output file
                                    outputFiles[idIndex] << kmerStr << " " << 1 << "\n";
                                    localKmersWritten++;
                                }
                            }

                            // Write genomes with multiple counts
                            for (const auto& countEntry : genomeInfo.genomesWithMultipleCounts) {
                                uint8_t idIndex = countEntry.first;
                                uint16_t count = countEntry.second;

                                if (idIndex >= startGenome && idIndex < endGenome) {
                                    // Write to the corresponding output file
                                    outputFiles[idIndex] << kmerStr << " " << count << "\n";
                                    localKmersWritten++;
                                }
                            }

                            localUniqueKmers++;
                        }
                    }

                    // Close all output files
                    for (auto& ofsEntry : outputFiles) {
                        ofsEntry.second.close();
                    }

                    // Update total kmers written
                    {
                        std::lock_guard<std::mutex> lock(countMutex);
                        totalKmersWritten += localKmersWritten;
                        uniqueKmers += localUniqueKmers;
                    }

                    std::cout << "[Writer " << threadId << "] Finished writing files. Kmers written: " << localKmersWritten << std::endl;
                });
            }

            // Wait for writer threads to finish
            for (auto& t : writerThreads) {
                t.join();
            }

            std::cout << "Total k-mers written to output files: " << totalKmersWritten << "\n";
        }

        auto end = std::chrono::high_resolution_clock::now();
        totalWritingTime += end - start;

        std::cout << "Total unique k-mers processed: " << uniqueKmers << "\n";
    }

    // Function to print timing information
    void printTimings() {
        std::cout << "\nProfiling Results:\n";
        std::cout << "Total sequences processed: " << totalSequencesProcessed << "\n";
        std::cout << "Total alignment blocks processed: " << totalAlignmentBlocks << "\n";
        std::cout << "Total k-mers processed (during processing): " << totalKMersProcessed << "\n";

        // Calculate total unique k-mers across partitions
        size_t totalUniqueKmers = 0;
        for (const auto& partition : kmerCountsPartitions) {
            totalUniqueKmers += partition.size();
        }
        std::cout << "Total unique k-mers: " << totalUniqueKmers << "\n";

        std::cout << "Time spent processing sequences: " << totalProcessSequenceTime.count() << " seconds\n";
        std::cout << "Time spent in merging KmerGroups: " << totalMergeTime.count() << " seconds\n";
        std::cout << "Time spent writing output files: " << totalWritingTime.count() << " seconds\n";
    }

private:
    int k;
    int numConsumers;
    bool useCanonicalKmers;
    bool producersDone;

    // Function to compute reverse complement
    KmerCodeType reverseComplement(KmerCodeType code, int k) {
        KmerCodeType revComp = 0;
        for (int i = 0; i < k; ++i) {
            uint8_t nucleotide = code & 0b11;  // Get the last two bits
            uint8_t complement = nucleotide ^ 0b11; // Complement the bits
            revComp = (revComp << 2) | complement;  // Append to revComp
            code >>= 2;  // Move to the next nucleotide
        }
        return revComp;
    }

    // Decode k-mer code to string remains unchanged
    std::string decodeKmer(KmerCodeType code) {
        std::string kmer(k, 'A');
        for (int i = k - 1; i >= 0; --i) {
            uint8_t bits = code & 0b11;
            switch (bits) {
                case 0b00: kmer[i] = 'A'; break;
                case 0b01: kmer[i] = 'C'; break;
                case 0b10: kmer[i] = 'G'; break;
                case 0b11: kmer[i] = 'T'; break;
            }
            code >>= 2;
        }
        return kmer;
    }
};

// Function to find chunk boundaries in the MAF file remains unchanged
void findMAFChunkBoundaries(const std::string& filename, int numChunks, std::vector<std::streampos>& chunkStartPositions) {
    std::ifstream mafFile(filename);
    if (!mafFile) {
        std::cerr << "Error opening MAF file: " << filename << std::endl;
        exit(1);
    }

    // Get file size
    mafFile.seekg(0, std::ios::end);
    std::streampos fileSize = mafFile.tellg();
    mafFile.seekg(0, std::ios::beg);

    // Estimate approximate chunk sizes
    std::vector<std::streampos> approxPositions(numChunks + 1);
    for (int i = 0; i <= numChunks; ++i) {
        approxPositions[i] = fileSize * i / numChunks;
    }

    // Adjust positions to start at the beginning of an alignment block ('a' line)
    chunkStartPositions.resize(numChunks + 1);
    chunkStartPositions[0] = 0; // First chunk starts at the beginning
    for (int i = 1; i < numChunks; ++i) {
        mafFile.seekg(approxPositions[i]);
        std::string line;
        while (std::getline(mafFile, line)) {
            if (line.empty()) continue;
            if (line[0] == 'a') {
                // Calculate the position at the start of the 'a' line
                std::streampos currentPos = mafFile.tellg();
                std::streamoff lineLength = static_cast<std::streamoff>(line.length()) + 1; // +1 for newline character
                chunkStartPositions[i] = currentPos - lineLength;
                break;
            }
            if (mafFile.tellg() >= fileSize) {
                chunkStartPositions[i] = fileSize;
                break;
            }
        }
    }
    chunkStartPositions[numChunks] = fileSize; // Last chunk ends at the end of file

    mafFile.close();
}

int main(int argc, char* argv[]) {
    auto programStartTime = std::chrono::high_resolution_clock::now();

    if (argc < 4) {
        std::cerr << "Usage: ./kmer_counter [options] <k-mer length> <MAF file> <number of threads>" << std::endl;
        std::cerr << "Options:\n"
                  << "  --complement, -c        Aggregate k-mers with their reverse complements\n"
                  << "  --single_file_output, -s   Write output to a single file\n";
        return 1;
    }

    // Variables for command-line arguments
    int k = 0;
    std::string mafFilename;
    int N = 0;
    bool useCanonicalKmers = false;
    bool singleFileOutput = false;

    // Parse command-line arguments
    std::vector<std::string> args(argv + 1, argv + argc);

    std::vector<std::string> positionalArgs;

    for (size_t i = 0; i < args.size(); ++i) {
        if (args[i] == "--complement" || args[i] == "-c") {
            useCanonicalKmers = true;
        } else if (args[i] == "--single_file_output" || args[i] == "-s") {
            singleFileOutput = true;
        } else {
            positionalArgs.push_back(args[i]);
        }
    }

    if (positionalArgs.size() != 3) {
        std::cerr << "Usage: ./kmer_counter [options] <k-mer length> <MAF file> <number of threads>" << std::endl;
        std::cerr << "Options:\n"
                  << "  --complement, -c        Aggregate k-mers with their reverse complements\n"
                  << "  --single_file_output, -s   Write output to a single file\n";
        return 1;
    }

    // Extract positional arguments
    k = std::stoi(positionalArgs[0]);
    mafFilename = positionalArgs[1];
    N = std::stoi(positionalArgs[2]);

    if (k <= 0 || k > 31) { // Limit k to 31 to fit in uint64_t (2 bits per base)
        std::cerr << "Error: k-mer length must be between 1 and 31." << std::endl;
        return 1;
    }

    if (N <= 1) {
        std::cerr << "Error: Number of threads must be greater than 1." << std::endl;
        return 1;
    }

    // Determine the number of producer and consumer threads
    int numProducers = N / 2;
    int numConsumers = N - numProducers;

    // Find chunk boundaries in the MAF file
    std::cout << "Finding chunk boundaries in the MAF file..." << std::endl;
    std::vector<std::streampos> chunkStartPositions;
    findMAFChunkBoundaries(mafFilename, numProducers, chunkStartPositions);

    // Choose KmerCodeType based on k-mer length
    if (k <= 16) {
        // Use uint32_t for k-mers up to length 16
        KmerCounter<uint32_t> counter(k, numConsumers, useCanonicalKmers);

        // Create per-consumer queues
        std::vector<moodycamel::ConcurrentQueue<KmerCounter<uint32_t>::KmerCountsMap>> kmerGroupQueues(numConsumers);

        // Create consumer threads
        std::vector<std::thread> consumerThreads;
        for (int i = 0; i < numConsumers; ++i) {
            consumerThreads.emplace_back(&KmerCounter<uint32_t>::consumeKmerGroups, &counter, std::ref(kmerGroupQueues[i]), i);
        }

        // Create producer threads
        std::vector<std::thread> producerThreads;
        for (int i = 0; i < numProducers; ++i) {
            std::streampos startPos = chunkStartPositions[i];
            std::streampos endPos = chunkStartPositions[i + 1];
            producerThreads.emplace_back(&KmerCounter<uint32_t>::processChunk, &counter, mafFilename, startPos, endPos, std::ref(kmerGroupQueues), i);
        }

        // Wait for producer threads to finish
        for (auto& t : producerThreads) {
            t.join();
        }

        // Signal the consumer threads that producers are done
        counter.setProducersDone();

        // Wait for consumer threads to finish
        for (auto& t : consumerThreads) {
            t.join();
        }

        // Write output files using N threads
        counter.writeOutputFiles(N, singleFileOutput);

        // Print timings
        counter.printTimings();
    } else {
        // Use uint64_t for k-mers up to length 31
        KmerCounter<uint64_t> counter(k, numConsumers, useCanonicalKmers);

        // Create per-consumer queues
        std::vector<moodycamel::ConcurrentQueue<KmerCounter<uint64_t>::KmerCountsMap>> kmerGroupQueues(numConsumers);

        // Create consumer threads
        std::vector<std::thread> consumerThreads;
        for (int i = 0; i < numConsumers; ++i) {
            consumerThreads.emplace_back(&KmerCounter<uint64_t>::consumeKmerGroups, &counter, std::ref(kmerGroupQueues[i]), i);
        }

        // Create producer threads
        std::vector<std::thread> producerThreads;
        for (int i = 0; i < numProducers; ++i) {
            std::streampos startPos = chunkStartPositions[i];
            std::streampos endPos = chunkStartPositions[i + 1];
            producerThreads.emplace_back(&KmerCounter<uint64_t>::processChunk, &counter, mafFilename, startPos, endPos, std::ref(kmerGroupQueues), i);
        }

        // Wait for producer threads to finish
        for (auto& t : producerThreads) {
            t.join();
        }

        // Signal the consumer threads that producers are done
        counter.setProducersDone();

        // Wait for consumer threads to finish
        for (auto& t : consumerThreads) {
            t.join();
        }

        // Write output files using N threads
        counter.writeOutputFiles(N, singleFileOutput);

        // Print timings
        counter.printTimings();
    }

    auto programEndTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> totalProgramTime = programEndTime - programStartTime;

    std::cout << "Total program execution time: " << totalProgramTime.count() << " seconds\n";

    return 0;
}
