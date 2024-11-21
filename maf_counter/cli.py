#!/usr/bin/env python3

import argparse
import subprocess
import sys
import os


def main():
    # Locate the maf_counter binary within the package
    script_dir = os.path.dirname(os.path.abspath(__file__))
    binary_path = os.path.join(script_dir, "bin", "maf_counter")

    if not os.path.isfile(binary_path):
        sys.stderr.write("Error: The 'maf_counter' binary is missing in the package.\n")
        sys.exit(1)

    # Define the argument parser
    parser = argparse.ArgumentParser(
        description="MAF Counter: A high-performance k-mer counting tool for MAF alignments."
    )
    parser.add_argument(
        "kmer_length",
        type=int,
        help="Length of the k-mers to count (integer between 1 and 31)."
    )
    parser.add_argument(
        "maf_file",
        type=str,
        help="Input MAF file."
    )
    parser.add_argument(
        "threads",
        type=int,
        help="Number of threads to use (integer > 1)."
    )
    parser.add_argument(
        "-c", "--complement",
        action="store_true",
        help="Aggregate k-mers with their reverse complements."
    )
    parser.add_argument(
        "-s", "--single_file_output",
        action="store_true",
        help="Write all k-mers to a single compressed file."
    )

    # Parse the arguments
    args = parser.parse_args()

    # Build the command for maf_counter
    command = [binary_path]

    if args.complement:
        command.append("-c")
    if args.single_file_output:
        command.append("-s")

    command.extend([str(args.kmer_length), args.maf_file, str(args.threads)])

    # Run the maf_counter binary
    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
        sys.stderr.write(f"Error: maf_counter exited with status {e.returncode}\n")
        sys.exit(e.returncode)


if __name__ == "__main__":
    main()
