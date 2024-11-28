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
#include <limits>      // For numeric_limits

#include "concurrentqueue.h" // For moodycamel::ConcurrentQueue
#include <boost/multiprecision/cpp_int.hpp>

using google::dense_hash_map;
using google::dense_hash_set;

// Define uint128_t using Boost.Multiprecision
namespace mp = boost::multiprecision;
using uint128_t = mp::uint128_t;
// Custom hash functions for uint16_t, uint32_t, uint64_t, and uint128_t keys
struct UInt16Hasher {
    std::size_t operator()(uint16_t key) const noexcept {
        return std::hash<uint16_t>{}(key);
    }
};
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
struct UInt128Hasher {
    std::size_t operator()(const uint128_t& key) const noexcept {
        // Combine the upper and lower 64 bits
        uint64_t lower = static_cast<uint64_t>(key & 0xFFFFFFFFFFFFFFFF);
        uint64_t upper = static_cast<uint64_t>(key >> 64);
        return std::hash<uint64_t>{}(lower) ^ std::hash<uint64_t>{}(upper);
    }
};
// Hash function selector based on KmerCodeType
template<typename T>
struct KmerCodeHasher;
template<>
struct KmerCodeHasher<uint16_t> : public UInt16Hasher {};
template<>
struct KmerCodeHasher<uint32_t> : public UInt32Hasher {};
template<>
struct KmerCodeHasher<uint64_t> : public UInt64Hasher {};
template<>
struct KmerCodeHasher<uint128_t> : public UInt128Hasher {};
// Enumeration for sequence types
enum class SequenceType { NUCLEOTIDES, AMINO_ACIDS };
// Custom struct to store genome counts in KmerGroups and KmerCounts
template<typename CountType, typename GenomeIdType>
struct GenomeInfo {
    dense_hash_set<GenomeIdType> genomesWithCountOne;
    dense_hash_map<GenomeIdType, CountType> genomesWithMultipleCounts;
    GenomeInfo() {
        genomesWithCountOne.set_empty_key(std::numeric_limits<GenomeIdType>::max());
        genomesWithCountOne.set_deleted_key(std::numeric_limits<GenomeIdType>::max() - 1);
        genomesWithMultipleCounts.set_empty_key(std::numeric_limits<GenomeIdType>::max() - 2);
        genomesWithMultipleCounts.set_deleted_key(std::numeric_limits<GenomeIdType>::max() - 3);
    }
    void increment(GenomeIdType genomeId) {
        if (genomesWithCountOne.find(genomeId) != genomesWithCountOne.end()) {
            // genomeId was already in genomesWithCountOne
            genomesWithCountOne.erase(genomeId);
            genomesWithMultipleCounts[genomeId] = static_cast<CountType>(2);
        } else if (genomesWithMultipleCounts.find(genomeId) != genomesWithMultipleCounts.end()) {
            // genomeId is in genomesWithMultipleCounts
            if (genomesWithMultipleCounts[genomeId] < std::numeric_limits<CountType>::max()) {
                genomesWithMultipleCounts[genomeId]++;
            }
            // else do nothing, count remains at max value
        } else {
            // genomeId is new
            genomesWithCountOne.insert(genomeId);
        }
    }
    // Function to get all genome IDs and counts
    std::vector<std::pair<GenomeIdType, CountType>> getGenomeIdCounts() const {
        std::vector<std::pair<GenomeIdType, CountType>> genomeIdCounts;
        genomeIdCounts.reserve(genomesWithCountOne.size() + genomesWithMultipleCounts.size());
        for (GenomeIdType genomeId : genomesWithCountOne) {
            genomeIdCounts.emplace_back(genomeId, static_cast<CountType>(1));
        }
        for (const auto& entry : genomesWithMultipleCounts) {
            genomeIdCounts.emplace_back(entry.first, entry.second);
        }
        return genomeIdCounts;
    }
};
template<typename KmerCodeType, typename CountType, typename GenomeIdType>
class KmerCounter {
public:
    // ID mappings
    std::unordered_map<std::string, GenomeIdType> idToIndex;
    std::vector<std::string> indexToId;
    std::mutex idMappingMutex; // Mutex for thread-safe genome ID mapping
    // K-mer counts
    typedef dense_hash_map<KmerCodeType, GenomeInfo<CountType, GenomeIdType>, KmerCodeHasher<KmerCodeType>> KmerCountsMap;
    std::vector<KmerCountsMap> kmerCountsPartitions;
    // Timing variables
    std::chrono::duration<double> totalProcessSequenceTime = std::chrono::duration<double>::zero();
    std::chrono::duration<double> totalMergeTime = std::chrono::duration<double>::zero();
    std::chrono::duration<double> totalWritingTime = std::chrono::duration<double>::zero();
    // Counters
    size_t totalKMersProcessed = 0;
    size_t totalSequencesProcessed = 0;
    size_t totalAlignmentBlocks = 0;
    KmerCounter(int k, int numConsumers, bool useCanonicalKmers, GenomeIdType maxGenomeId, SequenceType sequenceType)
        : k(k), numConsumers(numConsumers), useCanonicalKmers(useCanonicalKmers), maxGenomeId(maxGenomeId), sequenceType(sequenceType) {
        // Initialize kmerCountsPartitions
        kmerCountsPartitions.resize(numConsumers);
        for (auto& partition : kmerCountsPartitions) {
            partition.set_empty_key(std::numeric_limits<KmerCodeType>::max());
            partition.set_deleted_key(std::numeric_limits<KmerCodeType>::max() - 1);
        }
        producersDone = false;
        overflowDetected = false;
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
                        kg.set_empty_key(std::numeric_limits<KmerCodeType>::max());
                        kg.set_deleted_key(std::numeric_limits<KmerCodeType>::max() - 1);
                    }
                    size_t numSequences = sequences.size();
                    localTotalSequencesProcessed += numSequences;
                    // Process each sequence individually
                    for (const auto& seqPair : sequences) {
                        const std::string& id = seqPair.first;
                        const std::string& alignedSeq = seqPair.second;
                        GenomeIdType idIndex;
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
                                    if (idToIndex.size() >= maxGenomeId) {
                                        std::cerr << "Error: Exceeded maximum number of genome IDs (" << maxGenomeId << "). Re-run with --large_genome_count, -l option" << std::endl;
                                        exit(1);
                                    }
                                    idIndex = static_cast<GenomeIdType>(idToIndex.size());
                                    idToIndex[id] = idIndex;
                                    indexToId.push_back(id);
                                } else {
                                    idIndex = it->second;
                                }
                                // Lock is released when out of scope
                            }
                        }
                        if (sequenceType == SequenceType::NUCLEOTIDES) {
                            processNucleotideSequence(alignedSeq, idIndex, kmerGroups, localTotalKMersProcessed);
                        } else {
                            processAminoAcidSequence(alignedSeq, idIndex, kmerGroups, localTotalKMersProcessed);
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
                kg.set_empty_key(std::numeric_limits<KmerCodeType>::max());
                kg.set_deleted_key(std::numeric_limits<KmerCodeType>::max() - 1);
            }
            size_t numSequences = sequences.size();
            localTotalSequencesProcessed += numSequences;
            // Process each sequence individually
            for (const auto& seqPair : sequences) {
                const std::string& id = seqPair.first;
                const std::string& alignedSeq = seqPair.second;
                GenomeIdType idIndex;
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
                            if (idToIndex.size() >= maxGenomeId) {
                                std::cerr << "Error: Exceeded maximum number of genome IDs (" << maxGenomeId << ")." << std::endl;
                                exit(1);
                            }
                            idIndex = static_cast<GenomeIdType>(idToIndex.size());
                            idToIndex[id] = idIndex;
                            indexToId.push_back(id);
                        } else {
                            idIndex = it->second;
                        }
                        // Lock is released when out of scope
                    }
                }
                if (sequenceType == SequenceType::NUCLEOTIDES) {
                    processNucleotideSequence(alignedSeq, idIndex, kmerGroups, localTotalKMersProcessed);
                } else {
                    processAminoAcidSequence(alignedSeq, idIndex, kmerGroups, localTotalKMersProcessed);
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
                    const GenomeInfo<CountType, GenomeIdType>& blockGenomeInfo = entry.second;
                    auto& genomeInfo = kmerCountsPartitions[threadId][kmerCode];
                    // Merge genomesWithCountOne
                    for (GenomeIdType genomeId : blockGenomeInfo.genomesWithCountOne) {
                        if (genomeInfo.genomesWithCountOne.find(genomeId) != genomeInfo.genomesWithCountOne.end()) {
                            // Move to genomesWithMultipleCounts
                            genomeInfo.genomesWithCountOne.erase(genomeId);
                            genomeInfo.genomesWithMultipleCounts[genomeId] = static_cast<CountType>(2);
                        } else if (genomeInfo.genomesWithMultipleCounts.find(genomeId) != genomeInfo.genomesWithMultipleCounts.end()) {
                            // Increment count with overflow check
                            CountType& existingCount = genomeInfo.genomesWithMultipleCounts[genomeId];
                            if (existingCount < std::numeric_limits<CountType>::max()) {
                                existingCount += static_cast<CountType>(1);
                            }
                        } else {
                            genomeInfo.genomesWithCountOne.insert(genomeId);
                        }
                    }
                    // Merge genomesWithMultipleCounts
                    for (const auto& countEntry : blockGenomeInfo.genomesWithMultipleCounts) {
                        GenomeIdType genomeId = countEntry.first;
                        CountType count = countEntry.second;
                        if (genomeInfo.genomesWithCountOne.find(genomeId) != genomeInfo.genomesWithCountOne.end()) {
                            genomeInfo.genomesWithCountOne.erase(genomeId);
                            CountType newCount = static_cast<CountType>(1) + count;
                            if (newCount < count) { // Overflow check
                                newCount = std::numeric_limits<CountType>::max();
                                if (!overflowDetected){
                                    setOverFlowDetected();
                                }
                            }
                            genomeInfo.genomesWithMultipleCounts[genomeId] = newCount;
                        } else {
                            auto it = genomeInfo.genomesWithMultipleCounts.find(genomeId);
                            if (it != genomeInfo.genomesWithMultipleCounts.end()) {
                                CountType& existingCount = it->second;
                                if (existingCount > std::numeric_limits<CountType>::max() - count) {
                                    existingCount = std::numeric_limits<CountType>::max();
                                    if (!overflowDetected){
                                        setOverFlowDetected();
                                    }
                                } else {
                                    existingCount += count;
                                }
                            } else {
                                genomeInfo.genomesWithMultipleCounts[genomeId] = count;
                            }
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

    // Function to detect overflow and report it once
    void setOverFlowDetected(){
        overflowDetected = true;
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
                            const GenomeInfo<CountType, GenomeIdType>& genomeInfo = kmerEntry.second;

                            // Decode k-mer for output
                            std::string kmerStr = this->decodeKmer(kmerCode);

                            // Collect genome IDs and counts
                            std::vector<std::pair<GenomeIdType, CountType>> genomeIdCounts = genomeInfo.getGenomeIdCounts();

                            // Skip k-mers with no counts (shouldn't happen but added for safety)
                            if (genomeIdCounts.empty()) continue;

                            // Convert genome IDs to names and prepare output line
                            std::ostringstream oss;
                            oss << kmerStr;

                            for (const auto& idCountPair : genomeIdCounts) {
                                GenomeIdType idIndex = idCountPair.first;
                                CountType count = idCountPair.second;

                                // Ensure idIndex is within bounds
                                if (idIndex >= this->indexToId.size()) {
                                    std::cerr << "Error: Invalid idIndex " << static_cast<uint64_t>(idIndex) << std::endl;
                                    continue;
                                }

                                std::string genomeName = this->indexToId[idIndex];
                                oss << " " << genomeName << ": " << static_cast<uint32_t>(count);
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
            // Write to per-genome files
            size_t numGenomes = indexToId.size();
            size_t genomesPerThread = (numGenomes + numThreads - 1) / numThreads;

            std::vector<std::thread> writerThreads;

            for (int threadId = 0; threadId < numThreads; ++threadId) {
                size_t startGenome = threadId * genomesPerThread;
                size_t endGenome = std::min(startGenome + genomesPerThread, numGenomes);

                writerThreads.emplace_back([this, startGenome, endGenome, dirName, threadId, &totalKmersWritten, &uniqueKmers, &countMutex]() {
                    std::cout << "[Writer " << threadId << "] Started writing files for genomes " << startGenome << " to " << endGenome - 1 << std::endl;

                    // Open output files for assigned genomes
                    std::unordered_map<GenomeIdType, std::ofstream> outputFiles;
                    for (size_t i = startGenome; i < endGenome; ++i) {
                        std::string filename = std::string(dirName) + "/" + this->indexToId[i] + "_kmer_counts.txt";
                        outputFiles[static_cast<GenomeIdType>(i)].open(filename);
                        if (!outputFiles[static_cast<GenomeIdType>(i)]) {
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
                            const GenomeInfo<CountType, GenomeIdType>& genomeInfo = kmerEntry.second;

                            // Decode k-mer for output
                            std::string kmerStr = this->decodeKmer(kmerCode);

                            // Write genomes with count one
                            for (GenomeIdType idIndex : genomeInfo.genomesWithCountOne) {
                                if (idIndex >= startGenome && idIndex < endGenome) {
                                    // Write to the corresponding output file
                                    outputFiles[idIndex] << kmerStr << " " << static_cast<uint32_t>(1) << "\n";
                                    localKmersWritten++;
                                }
                            }

                            // Write genomes with multiple counts
                            for (const auto& countEntry : genomeInfo.genomesWithMultipleCounts) {
                                GenomeIdType idIndex = countEntry.first;
                                CountType count = countEntry.second;

                                if (idIndex >= startGenome && idIndex < endGenome) {
                                    // Write to the corresponding output file
                                    outputFiles[idIndex] << kmerStr << " " << static_cast<uint32_t>(count) << "\n";
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
        if (overflowDetected){
            std::cout << "WARNING: Max k-mer count reached. It is recommended to re-run with larger max_kmer_count (Options: 256, 65536, or 4294967296)\n";
        }
    }
private:
    int k;
    int numConsumers;
    bool useCanonicalKmers;
    bool producersDone;
    bool overflowDetected;
    GenomeIdType maxGenomeId;
    SequenceType sequenceType;
    // Function to compute reverse complement (only for nucleotides)
    KmerCodeType reverseComplement(KmerCodeType code, int k) {
        KmerCodeType revComp = 0;
        for (int i = 0; i < k; ++i) {
            uint8_t nucleotide = static_cast<uint8_t>(code & 0b11);  // Get the last two bits
            uint8_t complement = nucleotide ^ 0b11; // Complement the bits
            revComp = (revComp << 2) | complement;  // Append to revComp
            code >>= 2;  // Move to the next nucleotide
        }
        return revComp;
    }
    // Decode k-mer code to string
    std::string decodeKmer(KmerCodeType code) {
        if (sequenceType == SequenceType::NUCLEOTIDES) {
            std::string kmer(k, 'A');
            for (int i = k - 1; i >= 0; --i) {
                uint8_t bits = static_cast<uint8_t>(code & 0b11);
                switch (bits) {
                    case 0b00: kmer[i] = 'A'; break;
                    case 0b01: kmer[i] = 'C'; break;
                    case 0b10: kmer[i] = 'G'; break;
                    case 0b11: kmer[i] = 'T'; break;
                }
                code >>= 2;
            }
            return kmer;
        } else {
            // Amino acid decoding
            std::string kmer(k, 'A');
            for (int i = k - 1; i >= 0; --i) {
                uint8_t bits = static_cast<uint8_t>(code & 0b11111); // 5 bits
                char residue = codeToAminoAcid(bits);
                kmer[i] = residue;
                code >>= 5;
            }
            return kmer;
        }
    }
    
    // Map 5-bit codes back to amino acids
    char codeToAminoAcid(uint8_t code) {
        static const std::vector<char> codeToAminoAcidMap = {
            'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G',
            'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S',
            'T', 'W', 'Y', 'V', 'B', 'Z', 'U', 'X', 'O'
        };
        if (code < 25) {  // Changed from 20 to 25 to accommodate new codes
            return codeToAminoAcidMap[code];
        } else {
            return 'X'; // Unknown or invalid code
        }
    }

    // Mapping from amino acids to 5-bit codes
    const std::unordered_map<char, uint8_t> aminoAcidToCode = {
        {'A', 0}, {'R', 1}, {'N', 2}, {'D', 3}, {'C', 4},
        {'Q', 5}, {'E', 6}, {'G', 7}, {'H', 8}, {'I', 9},
        {'L', 10}, {'K', 11}, {'M', 12}, {'F', 13}, {'P', 14},
        {'S', 15}, {'T', 16}, {'W', 17}, {'Y', 18}, {'V', 19},
        {'B', 20}, {'Z', 21}, {'U', 22}, {'X', 23}, {'O', 24}  // Added new codes
    };
    // Function to process nucleotide sequences
    void processNucleotideSequence(const std::string& alignedSeq, GenomeIdType idIndex,
                                   std::vector<KmerCountsMap>& kmerGroups, size_t& localTotalKMersProcessed) {
        size_t seqLength = alignedSeq.size();
        KmerCodeType kmerCode = 0;
        size_t validBases = 0;
        KmerCodeType mask = (k >= static_cast<int>(sizeof(KmerCodeType) * 4)) ? static_cast<KmerCodeType>(-1)
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
                int partitionIndex = static_cast<int>(finalKmerCode % numConsumers);
                auto& genomeInfo = kmerGroups[partitionIndex][finalKmerCode];
                genomeInfo.increment(idIndex);
                localTotalKMersProcessed++;
            }
        }
    }
    // Function to process amino acid sequences
    void processAminoAcidSequence(const std::string& alignedSeq, GenomeIdType idIndex,
                                  std::vector<KmerCountsMap>& kmerGroups, size_t& localTotalKMersProcessed) {
        size_t seqLength = alignedSeq.size();
        KmerCodeType kmerCode = 0;
        size_t validResidues = 0;
        KmerCodeType mask = ((static_cast<KmerCodeType>(1) << (5 * k)) - 1);
        for (size_t i = 0; i < seqLength; ++i) {
            char residue = alignedSeq[i];
            residue = toupper(residue);
            if (residue == '-') {
                continue; // Skip gaps
            }
            auto it = aminoAcidToCode.find(residue);
            if (it == aminoAcidToCode.end()) {
                // Invalid amino acid, reset k-mer
                validResidues = 0;
                kmerCode = 0;
                continue;
            }
            uint8_t code = it->second;
            kmerCode = ((kmerCode << 5) | code) & mask;
            validResidues++;
            if (validResidues >= static_cast<size_t>(k)) {
                KmerCodeType finalKmerCode = kmerCode;
                // No reverse complement for amino acids
                int partitionIndex = static_cast<int>(finalKmerCode % numConsumers);
                auto& genomeInfo = kmerGroups[partitionIndex][finalKmerCode];
                genomeInfo.increment(idIndex);
                localTotalKMersProcessed++;
            }
        }
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
    
    // Template function to run the KmerCounter with specified types
    template<typename KmerCodeType, typename CountType, typename GenomeIdType>
    void runKmerCounter(int k, int numConsumers, bool useCanonicalKmers, const std::string& mafFilename, int numProducers,
                        const std::vector<std::streampos>& chunkStartPositions, bool singleFileOutput, GenomeIdType maxGenomeId,
                        SequenceType sequenceType) {
        KmerCounter<KmerCodeType, CountType, GenomeIdType> counter(k, numConsumers, useCanonicalKmers, maxGenomeId, sequenceType);
    
        // Create per-consumer queues
        std::vector<moodycamel::ConcurrentQueue<typename KmerCounter<KmerCodeType, CountType, GenomeIdType>::KmerCountsMap>> kmerGroupQueues(numConsumers);
    
        // Create consumer threads
        std::vector<std::thread> consumerThreads;
        for (int i = 0; i < numConsumers; ++i) {
            consumerThreads.emplace_back(&KmerCounter<KmerCodeType, CountType, GenomeIdType>::consumeKmerGroups, &counter, std::ref(kmerGroupQueues[i]), i);
        }
    
        // Create producer threads
        std::vector<std::thread> producerThreads;
        for (int i = 0; i < numProducers; ++i) {
            std::streampos startPos = chunkStartPositions[i];
            std::streampos endPos = chunkStartPositions[i + 1];
            producerThreads.emplace_back(&KmerCounter<KmerCodeType, CountType, GenomeIdType>::processChunk, &counter, mafFilename, startPos, endPos, std::ref(kmerGroupQueues), i);
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
        counter.writeOutputFiles(numConsumers + numProducers, singleFileOutput);
    
        // Print timings
        counter.printTimings();
    }
    
    // Enums for CountType and GenomeIdType
    enum class CountTypeEnum { UINT8, UINT16, UINT32 };
    enum class GenomeIdTypeEnum { UINT8, UINT16 };
    enum class KmerCodeTypeEnum { UINT16, UINT32, UINT64, UINT128 };
    
    // Function to determine KmerCodeTypeEnum based on k-mer length and sequence type
    KmerCodeTypeEnum determineKmerCodeType(int k, SequenceType sequenceType) {
        if (sequenceType == SequenceType::NUCLEOTIDES) {
            if (k <= 16) {
                return KmerCodeTypeEnum::UINT32;
            } else if (k <= 32) {
                return KmerCodeTypeEnum::UINT64;
            } else {
                return KmerCodeTypeEnum::UINT128;
            }
        } else { // AMINO_ACIDS
            if (k == 3) {
                return KmerCodeTypeEnum::UINT16;
            } else if (k >= 4 && k <= 6) {
                return KmerCodeTypeEnum::UINT32;
            } else if (k >= 7 && k <= 12) {
                return KmerCodeTypeEnum::UINT64;
            } else if (k >= 13 && k <= 25) {
                return KmerCodeTypeEnum::UINT128;
            } else {
                std::cerr << "Error: Invalid k-mer length for amino acids. Supported k values are 3 to 25." << std::endl;
                exit(1);
            }
        }
    }
    int main(int argc, char* argv[]) {
    auto programStartTime = std::chrono::high_resolution_clock::now();
    if (argc < 4) {
        std::cerr << "Usage: ./kmer_counter [options] <k-mer length> <MAF file> <number of threads>" << std::endl;
        std::cerr << "Options:\n"
                  << "  --complement, -c             Aggregate k-mers with their reverse complements\n"
                  << "  --single_file_output, -s     Write output to a single file\n"
                  << "  --max_kmer_count <value>, -m <value>   Set maximum k-mer count (256, 65536, or 4294967296)\n"
                  << "  --large_genome_count, -l     Support more than 256 genomes (up to 65,536)\n"
                  << "  --sequence_type <type>, -t <type>      Set sequence type ('nucleotides' or 'amino_acids')\n";
        return 1;
    }
    // Variables for command-line arguments
    int k = 0;
    std::string mafFilename;
    int N = 0;
    bool useCanonicalKmers = false;
    bool singleFileOutput = false;
    uint64_t maxKmerCount = 65536; // Default value
    bool largeGenomeCount = false;
    std::string sequenceTypeStr = "nucleotides"; // Default sequence type
    // Parse command-line arguments
    std::vector<std::string> args(argv + 1, argv + argc);
    std::vector<std::string> positionalArgs;
    for (size_t i = 0; i < args.size(); ++i) {
        if (args[i] == "--complement" || args[i] == "-c") {
            useCanonicalKmers = true;
        } else if (args[i] == "--single_file_output" || args[i] == "-s") {
            singleFileOutput = true;
        } else if (args[i] == "--max_kmer_count" || args[i] == "-m") {
            if (i + 1 >= args.size()) {
                std::cerr << "Error: --max_kmer_count option requires an argument." << std::endl;
                return 1;
            }
            maxKmerCount = std::stoull(args[++i]);
        } else if (args[i] == "--large_genome_count" || args[i] == "-l") {
            largeGenomeCount = true;
        } else if (args[i] == "--sequence_type" || args[i] == "-t") {
            if (i + 1 >= args.size()) {
                std::cerr << "Error: --sequence_type option requires an argument ('nucleotides' or 'amino_acids')." << std::endl;
                return 1;
            }
            sequenceTypeStr = args[++i];
            if (sequenceTypeStr != "nucleotides" && sequenceTypeStr != "amino_acids") {
                std::cerr << "Error: Invalid sequence type. Must be 'nucleotides' or 'amino_acids'." << std::endl;
                return 1;
            }
        } else {
            positionalArgs.push_back(args[i]);
        }
    }
    if (positionalArgs.size() != 3) {
        std::cerr << "Usage: ./kmer_counter [options] <k-mer length> <MAF file> <number of threads>" << std::endl;
        std::cerr << "Options:\n"
                  << "  --complement, -c             Aggregate k-mers with their reverse complements\n"
                  << "  --single_file_output, -s     Write output to a single file\n"
                  << "  --max_kmer_count <value>, -m <value>   Set maximum k-mer count (256, 65536, or 4294967296)\n"
                  << "  --large_genome_count, -l     Support more than 256 genomes (up to 65,536)\n"
                  << "  --sequence_type <type>, -t <type>      Set sequence type ('nucleotides' or 'amino_acids')\n";
        return 1;
    }
    // Extract positional arguments
    k = std::stoi(positionalArgs[0]);
    mafFilename = positionalArgs[1];
    N = std::stoi(positionalArgs[2]);
    // Determine sequence type
    SequenceType sequenceType = (sequenceTypeStr == "nucleotides") ? SequenceType::NUCLEOTIDES : SequenceType::AMINO_ACIDS;
    // Validate k-mer length
    if (sequenceType == SequenceType::NUCLEOTIDES) {
        if (k <= 0 || k > 63) { // Limit k to 63 to fit in uint128_t (2 bits per base)
            std::cerr << "Error: k-mer length must be between 1 and 63 for nucleotides." << std::endl;
            return 1;
        }
    } else {
        if (k < 3 || k > 25) {
            std::cerr << "Error: k-mer length must be between 3 and 25 for amino acids." << std::endl;
            return 1;
        }
        if (useCanonicalKmers) {
            std::cerr << "Error: Canonical Kmers option incompatible with --sequence_type, -t amino_acids." << std::endl;
            return 1;
        }
    }
    if (N <= 1) {
        std::cerr << "Error: Number of threads must be greater than 1." << std::endl;
        return 1;
    }
    if (maxKmerCount != 256 && maxKmerCount != 65536 && maxKmerCount != 4294967296) {
        std::cerr << "Error: max_kmer_count must be one of 256, 65536, or 4294967296." << std::endl;
        return 1;
    }
    // Determine the CountType based on maxKmerCount
    CountTypeEnum countType;
    if (maxKmerCount <= 256) {
        countType = CountTypeEnum::UINT8;
    } else if (maxKmerCount <= 65536) {
        countType = CountTypeEnum::UINT16;
    } else {
        countType = CountTypeEnum::UINT32;
    }
    // Determine GenomeIdType
    GenomeIdTypeEnum genomeIdType = largeGenomeCount ? GenomeIdTypeEnum::UINT16 : GenomeIdTypeEnum::UINT8;
    // Determine the number of producer and consumer threads
    int numProducers = N / 2;
    int numConsumers = N - numProducers;
    // Find chunk boundaries in the MAF file
    std::cout << "Finding chunk boundaries in the MAF file..." << std::endl;
    std::vector<std::streampos> chunkStartPositions;
    findMAFChunkBoundaries(mafFilename, numProducers, chunkStartPositions);
    // Determine KmerCodeTypeEnum
    KmerCodeTypeEnum kmerCodeTypeEnum = determineKmerCodeType(k, sequenceType);
    // Choose KmerCodeType based on k-mer length and sequence type
    if (kmerCodeTypeEnum == KmerCodeTypeEnum::UINT16) {
        using KmerCodeType = uint16_t;
        if (countType == CountTypeEnum::UINT8) {
            if (genomeIdType == GenomeIdTypeEnum::UINT8) {
                runKmerCounter<KmerCodeType, uint8_t, uint8_t>(
                    k, numConsumers, useCanonicalKmers, mafFilename, numProducers,
                    chunkStartPositions, singleFileOutput, static_cast<uint8_t>(254), sequenceType);
            } else {
                runKmerCounter<KmerCodeType, uint8_t, uint16_t>(
                    k, numConsumers, useCanonicalKmers, mafFilename, numProducers,
                    chunkStartPositions, singleFileOutput, static_cast<uint16_t>(65534), sequenceType);
            }
        } else if (countType == CountTypeEnum::UINT16) {
            if (genomeIdType == GenomeIdTypeEnum::UINT8) {
                runKmerCounter<KmerCodeType, uint16_t, uint8_t>(
                    k, numConsumers, useCanonicalKmers, mafFilename, numProducers,
                    chunkStartPositions, singleFileOutput, static_cast<uint8_t>(254), sequenceType);
            } else {
                runKmerCounter<KmerCodeType, uint16_t, uint16_t>(
                    k, numConsumers, useCanonicalKmers, mafFilename, numProducers,
                    chunkStartPositions, singleFileOutput, static_cast<uint16_t>(65534), sequenceType);
            }
        } else { // UINT32
            if (genomeIdType == GenomeIdTypeEnum::UINT8) {
                runKmerCounter<KmerCodeType, uint32_t, uint8_t>(
                    k, numConsumers, useCanonicalKmers, mafFilename, numProducers,
                    chunkStartPositions, singleFileOutput, static_cast<uint8_t>(254), sequenceType);
            } else {
                runKmerCounter<KmerCodeType, uint32_t, uint16_t>(
                    k, numConsumers, useCanonicalKmers, mafFilename, numProducers,
                    chunkStartPositions, singleFileOutput, static_cast<uint16_t>(65534), sequenceType);
            }
        }
    } else if (kmerCodeTypeEnum == KmerCodeTypeEnum::UINT32) {
        using KmerCodeType = uint32_t;
        if (countType == CountTypeEnum::UINT8) {
            if (genomeIdType == GenomeIdTypeEnum::UINT8) {
                runKmerCounter<KmerCodeType, uint8_t, uint8_t>(
                    k, numConsumers, useCanonicalKmers, mafFilename, numProducers,
                    chunkStartPositions, singleFileOutput, static_cast<uint8_t>(254), sequenceType);
            } else {
                runKmerCounter<KmerCodeType, uint8_t, uint16_t>(
                    k, numConsumers, useCanonicalKmers, mafFilename, numProducers,
                    chunkStartPositions, singleFileOutput, static_cast<uint16_t>(65534), sequenceType);
            }
        } else if (countType == CountTypeEnum::UINT16) {
            if (genomeIdType == GenomeIdTypeEnum::UINT8) {
                runKmerCounter<KmerCodeType, uint16_t, uint8_t>(
                    k, numConsumers, useCanonicalKmers, mafFilename, numProducers,
                    chunkStartPositions, singleFileOutput, static_cast<uint8_t>(254), sequenceType);
            } else {
                runKmerCounter<KmerCodeType, uint16_t, uint16_t>(
                    k, numConsumers, useCanonicalKmers, mafFilename, numProducers,
                    chunkStartPositions, singleFileOutput, static_cast<uint16_t>(65534), sequenceType);
            }
        } else { // UINT32
            if (genomeIdType == GenomeIdTypeEnum::UINT8) {
                runKmerCounter<KmerCodeType, uint32_t, uint8_t>(
                    k, numConsumers, useCanonicalKmers, mafFilename, numProducers,
                    chunkStartPositions, singleFileOutput, static_cast<uint8_t>(254), sequenceType);
            } else {
                runKmerCounter<KmerCodeType, uint32_t, uint16_t>(
                    k, numConsumers, useCanonicalKmers, mafFilename, numProducers,
                    chunkStartPositions, singleFileOutput, static_cast<uint16_t>(65534), sequenceType);
            }
        }
    } else if (kmerCodeTypeEnum == KmerCodeTypeEnum::UINT64) {
        using KmerCodeType = uint64_t;
        if (countType == CountTypeEnum::UINT8) {
            if (genomeIdType == GenomeIdTypeEnum::UINT8) {
                runKmerCounter<KmerCodeType, uint8_t, uint8_t>(
                    k, numConsumers, useCanonicalKmers, mafFilename, numProducers,
                    chunkStartPositions, singleFileOutput, static_cast<uint8_t>(254), sequenceType);
            } else {
                runKmerCounter<KmerCodeType, uint8_t, uint16_t>(
                    k, numConsumers, useCanonicalKmers, mafFilename, numProducers,
                    chunkStartPositions, singleFileOutput, static_cast<uint16_t>(65534), sequenceType);
            }
        } else if (countType == CountTypeEnum::UINT16) {
            if (genomeIdType == GenomeIdTypeEnum::UINT8) {
                runKmerCounter<KmerCodeType, uint16_t, uint8_t>(
                    k, numConsumers, useCanonicalKmers, mafFilename, numProducers,
                    chunkStartPositions, singleFileOutput, static_cast<uint8_t>(254), sequenceType);
            } else {
                runKmerCounter<KmerCodeType, uint16_t, uint16_t>(
                    k, numConsumers, useCanonicalKmers, mafFilename, numProducers,
                    chunkStartPositions, singleFileOutput, static_cast<uint16_t>(65534), sequenceType);
            }
        } else { // UINT32
            if (genomeIdType == GenomeIdTypeEnum::UINT8) {
                runKmerCounter<KmerCodeType, uint32_t, uint8_t>(
                    k, numConsumers, useCanonicalKmers, mafFilename, numProducers,
                    chunkStartPositions, singleFileOutput, static_cast<uint8_t>(254), sequenceType);
            } else {
                runKmerCounter<KmerCodeType, uint32_t, uint16_t>(
                    k, numConsumers, useCanonicalKmers, mafFilename, numProducers,
                    chunkStartPositions, singleFileOutput, static_cast<uint16_t>(65534), sequenceType);
            }
        }
    } else if (kmerCodeTypeEnum == KmerCodeTypeEnum::UINT128) {
        using KmerCodeType = uint128_t;
        if (countType == CountTypeEnum::UINT8) {
            if (genomeIdType == GenomeIdTypeEnum::UINT8) {
                runKmerCounter<KmerCodeType, uint8_t, uint8_t>(
                    k, numConsumers, useCanonicalKmers, mafFilename, numProducers,
                    chunkStartPositions, singleFileOutput, static_cast<uint8_t>(254), sequenceType);
            } else {
                runKmerCounter<KmerCodeType, uint8_t, uint16_t>(
                    k, numConsumers, useCanonicalKmers, mafFilename, numProducers,
                    chunkStartPositions, singleFileOutput, static_cast<uint16_t>(65534), sequenceType);
            }
        } else if (countType == CountTypeEnum::UINT16) {
            if (genomeIdType == GenomeIdTypeEnum::UINT8) {
                runKmerCounter<KmerCodeType, uint16_t, uint8_t>(
                    k, numConsumers, useCanonicalKmers, mafFilename, numProducers,
                    chunkStartPositions, singleFileOutput, static_cast<uint8_t>(254), sequenceType);
            } else {
                runKmerCounter<KmerCodeType, uint16_t, uint16_t>(
                    k, numConsumers, useCanonicalKmers, mafFilename, numProducers,
                    chunkStartPositions, singleFileOutput, static_cast<uint16_t>(65534), sequenceType);
            }
        } else { // UINT32
            if (genomeIdType == GenomeIdTypeEnum::UINT8) {
                runKmerCounter<KmerCodeType, uint32_t, uint8_t>(
                    k, numConsumers, useCanonicalKmers, mafFilename, numProducers,
                    chunkStartPositions, singleFileOutput, static_cast<uint8_t>(254), sequenceType);
            } else {
                runKmerCounter<KmerCodeType, uint32_t, uint16_t>(
                    k, numConsumers, useCanonicalKmers, mafFilename, numProducers,
                    chunkStartPositions, singleFileOutput, static_cast<uint16_t>(65534), sequenceType);
            }
        }
    } else {
        std::cerr << "Error: Unsupported KmerCodeType." << std::endl;
        return 1;
    }
    auto programEndTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> totalProgramTime = programEndTime - programStartTime;
    std::cout << "Total program execution time: " << totalProgramTime.count() << " seconds\n";
    return 0;
}
    
