#!/usr/bin/env python
"""
samplesheet.py -m manifest.csv -f folder1,folder2
read all paired fastq files in folders and create a samplesheet.csv file.
"""

import argparse
import csv
import os
from collections import defaultdict


def get_manifest(manifest_file):
    """
    manifest_file
    ```
    sample,prefix
    sample_X,prefix_X
    sample_Y,prefix_Y
    ```
    returns:
        A dict of {prefix: sample}
    """
    manifest = {}
    with open(manifest_file) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            manifest[row["prefix"]] = row["sample"]
    return manifest


def get_pair(file_name):
    # determine if it is R1 or R2
    for x in (1, 2):
        for y in (f"_{x}_", f"_{x}.", f"_R{x}_", f"_R{x}."):
            if y in file_name:
                return x


def find_fastq_files(folders, manifest):
    fastq_files = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    for folder in folders:
        folder = os.path.abspath(folder)
        for root, _, files in os.walk(folder):
            for file in files:
                if not file.endswith(".gz"):
                    continue
                for prefix, sample in manifest.items():
                    if prefix + "_" in file:
                        x = get_pair(file)
                        # absolute path
                        path = os.path.join(root, file)
                        fastq_files[sample][prefix][x].append(path)
        print(f"Found {len(fastq_files)} samples in folder {folder}")
    return fastq_files


def write_samplesheet(samplesheet_file, fastq_files, manifest):
    n = 0
    with open(samplesheet_file, "w", newline="") as csvfile:
        fieldnames = ["sample", "fastq_1", "fastq_2"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for prefix, sample in manifest.items():
            r1, r2 = fastq_files[sample][prefix][1], fastq_files[sample][prefix][2]
            if len(r1) != len(r2):
                raise ValueError(
                    f"Number of R1 files ({len(r1)}) and R2 files ({len(r2)}) do not match for sample {sample} and prefix {prefix}"
                )
            if not r1:
                print(f"WARNING: No files found for sample {sample} and prefix {prefix}")
                continue
            for f1, f2 in zip(sorted(r1), sorted(r2)):
                writer.writerow({"sample": sample, "fastq_1": f1, "fastq_2": f2})
            n += 1
    return n


def main():
    parser = argparse.ArgumentParser(description="Generate samplesheet from fastq files.")
    parser.add_argument(
        "-m", "--manifest", required=True, help="Path to the manifest CSV file containing prefix-sample mapping."
    )
    parser.add_argument(
        "-f", "--folders", required=True, help="Comma-separated paths to folders to search for fastq files."
    )
    parser.add_argument("-o", "--output", help="Output samplesheet file.", default="samplesheet.csv")

    args = parser.parse_args()

    folders = args.folders.split(",")
    manifest = get_manifest(args.manifest)
    print(f"Found {len(manifest)} samples in {args.manifest}")
    fastq_files = find_fastq_files(folders, manifest)
    n = write_samplesheet(args.output, fastq_files, manifest)
    print(f"Wrote {n} samples to {args.output}")


if __name__ == "__main__":
    main()