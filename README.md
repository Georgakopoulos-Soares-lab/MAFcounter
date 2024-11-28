# MAF Counter

MAF Counter is a multithreaded tool designed to efficiently extract and count k-mers from multiple genome alignments in MAF format. It utilizes producer-consumer threading with concurrent queues to parallelize k-mer extraction and aggregation. The tool supports canonical/reverse-complement k-mer handling and provides flexible output options, including per-genome k-mer files and a consolidated single-file format.


## Algorithm Overview
1. **Chunk Division**: The input MAF file is divided into approximately equal-sized chunks, each containing whole alignment blocks, ensuring no block spans multiple chunks.
2. **Parallel Processing**: Producer threads process each chunk independently, using a sliding window to extract k-mers and their counts.
3. **Consumer Merging**: Consumer threads merge k-mer groups into a final aggregated data structure, avoiding conflicts by partitioning based on hashed k-mer keys.
4. **File Writing**: Once processing is complete, output is written either as per-genome files in Jellyfish format or a single compressed file if the `-s` flag is used.

## Install MAFcounter CLI python wrapper through pip
```bash
pip install maf-counter
maf_counter --help
```

## Download and Execute MAFcounter binary from releases
To 
```bash
wget https://github.com/Georgakopoulos-Soares-lab/MAFcounter/releases/download/v1.0_linux_x86/maf_counter
./maf_counter 15 input.maf 16
```

## Compilation
To compile MAF Counter, use the following command after downloading 
1) Google Sparse Hash ( https://github.com/sparsehash/sparsehash )
2) Concurrent Queue ( https://github.com/cameron314/concurrentqueue )
3) Boost Multiprecision CPPInt https://www.boost.org/doc/libs/1_86_0/libs/multiprecision/doc/html/index.html
```bash
g++ -std=c++11 -O3 -o maf_counter maf_counter.cpp -I /path/to/google_sparse_hash -I /path/to/concurrent_queue -I /path/to/boost -pthread -lrt
```

## Usage
```
./maf_counter [options] <k-mer length> <MAF file> <number of threads>
```

## Options
```
-c, --complement: Aggregate k-mers with their reverse complements.
-s, --single_file_output: Write all k-mers to a single compressed file.
-m, --max_kmer_count: Set maximum k-mer count (accepts values: 256, 65536, or 4294967296).
-l, --large_genome_count: Support more than 256 genomes (up to 65,536).
-t, --sequence_type: Set sequence type (accepts values: 'nucleotides' or 'amino_acids').

```

## Examples

- Extract 15mers from input.maf using 8 producers and 8 consumers ( suitable ideally for 16 processor cores )
```
./maf_counter 15 input.maf 16
```
- The same options but aggregate each kmer with its reverse complement (writing the lexicographically first) and output the results in sinle file mode
```
./maf_counter -c -s 15 input.maf 16
```
- Extract 15mers from maf file with more than 255 genomes using 10 cores.
```
./maf_counter 15 more_genomes.maf 10 --large_genome_count
```
- Extract 5mers from proteome instead of genome while keeping the max kmer count to 256 for lower memory usage.
```
./maf_counter 5 example_proteome.maf 5 --max_kmer_count 256 -t amino_acids
```

## Example Files
Example files can be found in the [maf_counter/example_files](https://github.com/Georgakopoulos-Soares-lab/MAFcounter/tree/master/maf_counter/example_files) directory and they include a synthetic genome and a synthetic proteome for tool demonstration



## Output Format
By default, the tool generates per-genome k-mer files in the results_counter directory. Each file is in Jellyfish format, containing k-mers and their counts for the corresponding genome ID.

Example file (<b>genome1_kmer_counts.txt</b>):
```yaml
ATCGG  1401
TTGGC  1233
```
If the -s option is used, a single output file is created in the format:
```yaml
ATCGG genomeid1:1401, genomeid2:1200
TTGGC genomeid1:1233, genomeid3:4123
```

## License
This project is licensed under the [GNU GPL v3](LICENSE).

## Contact

For any questions or support, please contact 
- izg5139@psu.edu
- mpp5977@psu.edu 
- kap6605@psu.edu
- ioannis.mouratidis@psu.edu
