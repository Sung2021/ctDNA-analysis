#!/usr/bin/env python3
"""
02_alignment.py
BWA-MEM2 alignment → sorted BAM generation.

ctDNA-specific assumption:
    - Input FASTQ files are already UMI-extracted (UMI-tools extract or equivalent).
    - Therefore each sample has paired reads: R1 (sequence), R2 (sequence or UMI-carried read depending on pipeline).
    - No UMI processing is performed here; pure alignment only.

Usage:
    python 02_alignment.py \
        --fastq_dir ./data/processed/umi_extracted \
        --reference ./refs/hg38.fa \
        --out_dir ./results/alignment \
        --threads 8
"""

import os
import argparse
import glob
import subprocess

def run_cmd(cmd):
    print("[CMD] " + " ".join(cmd))
    subprocess.run(cmd, check=True)

def find_fastq_pairs(fastq_dir):
    """
    Find paired FASTQ files assuming naming *_R1*.fastq.gz and *_R2*.fastq.gz.
    This code assumes UMI extraction is already completed.
    """
    r1_files = sorted(glob.glob(os.path.join(fastq_dir, "*_R1*.fastq.gz")))
    pairs = []

    for r1 in r1_files:
        r2 = r1.replace("_R1", "_R2")
        if os.path.exists(r2):
            sample = os.path.basename(r1).split("_R1")[0]
            pairs.append((sample, r1, r2))
        else:
            print(f"[WARN] R2 not found for {r1}")

    if len(pairs) == 0:
        raise RuntimeError("No FASTQ pairs found.")
    return pairs

def run_alignment(pairs, reference, out_dir, threads):
    os.makedirs(out_dir, exist_ok=True)

    for sample, r1, r2 in pairs:
        bam = os.path.join(out_dir, f"{sample}.sorted.bam")

        # BWA-MEM2 alignment command
        cmd_align = [
            "bwa-mem2", "mem",
            "-t", str(threads),
            reference, r1, r2
        ]

        # Pipe → samtools sort
        cmd_samtools_sort = [
            "samtools", "sort",
            "-@", str(threads),
            "-o", bam,
            "-"
        ]

        print(f"[INFO] Processing sample: {sample}")

        p1 = subprocess.Popen(cmd_align, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(cmd_samtools_sort, stdin=p1.stdout)
        p2.communicate()
        p1.wait()

        if p2.returncode != 0:
            raise RuntimeError(f"Sorting failed for {sample}")

        run_cmd(["samtools", "index", bam])
        print(f"[INFO] Completed: {bam}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fastq_dir", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()

    pairs = find_fastq_pairs(args.fastq_dir)
    run_alignment(pairs, args.reference, args.out_dir, args.threads)

if __name__ == "__main__":
    main()
