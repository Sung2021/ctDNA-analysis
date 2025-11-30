#!/usr/bin/env python3
"""
03_umi_processing_umitools.py
UMI-tools deduplication for ctDNA.

Assumptions:
    - Input BAM files come from 02_alignment.py (sorted + indexed + RG added).
    - UMI extraction was already done upstream (FASTQ → UMI-tools extract).
    - UMIs are stored in RX tag (standard UMI-tools extract behavior).

Usage:
    python 03_umi_processing_umitools.py \
        --bam_dir ./results/alignment \
        --out_dir ./results/umi_dedup \
        --method directional \
        --threads 4
"""

import os
import argparse
import glob
import subprocess

def run_cmd(cmd):
    print("[CMD] " + " ".join(cmd))
    subprocess.run(cmd, check=True)

def find_bam_files(bam_dir):
    bams = sorted(glob.glob(os.path.join(bam_dir, "*.sorted.bam")))
    if len(bams) == 0:
        raise FileNotFoundError(f"No sorted BAM files found in: {bam_dir}")
    return bams

def run_umi_dedup(bams, out_dir, method, threads):
    os.makedirs(out_dir, exist_ok=True)

    for bam in bams:
        sample = os.path.basename(bam).replace(".sorted.bam", "")
        out_bam = os.path.join(out_dir, f"{sample}.dedup.bam")
        metrics = os.path.join(out_dir, f"{sample}.dedup_metrics.txt")

        # UMI-tools dedup command
        cmd = [
            "umi_tools", "dedup",
            "--extract-umi-method=tag",
            "--umi-tag=RX",
            "--method", method,
            "--threads", str(threads),
            "--stdin", bam,
            "--stdout", out_bam,
            "--log", metrics,
        ]

        print(f"[INFO] Deduplicating sample: {sample}")
        run_cmd(cmd)

        # Index BAM with threads
        run_cmd(["samtools", "index", "-@", str(threads), out_bam])
        print(f"[INFO] Completed: {out_bam}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam_dir", required=True, help="Input directory containing *.sorted.bam")
    parser.add_argument("--out_dir", required=True, help="Output directory for deduplicated BAM files")
    parser.add_argument("--method", default="directional",
                        choices=["unique", "percentile", "directional"],
                        help="UMI-tools dedup method")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads to use")
    args = parser.parse_args()

    # Step 1: find aligned BAMs
    bam_files = find_bam_files(args.bam_dir)

    # Step 2: perform UMI dedup
    run_umi_dedup(bam_files, args.out_dir, args.method, args.threads)

if __name__ == "__main__":
    main()
