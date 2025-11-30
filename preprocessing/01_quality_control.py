#!/usr/bin/env python3
"""
01_quality_control.py
FASTQ QC using FastQC + summary report generation.

Usage:
    python 01_quality_control.py --fastq_dir ./data/raw --out_dir ./results/qc
"""

import os
import argparse
import subprocess
import glob
import pandas as pd

def run_fastqc(fastq_files, out_dir, threads=4):
    os.makedirs(out_dir, exist_ok=True)
    for fq in fastq_files:
        cmd = [
            "fastqc",
            fq,
            "-t", str(threads),
            "-o", out_dir
        ]
        print(f"[INFO] Running FastQC on {os.path.basename(fq)}")
        subprocess.run(cmd, check=True)

def parse_fastqc_summary(fastqc_dir):
    """
    Collects all 'summary.txt' files from FastQC outputs
    and builds a merged QC table.
    """
    records = []

    summary_files = glob.glob(os.path.join(fastqc_dir, "*_fastqc", "summary.txt"))
    for sf in summary_files:
        sample = os.path.basename(os.path.dirname(sf)).replace("_fastqc", "")
        with open(sf) as f:
            for line in f:
                status, metric, _ = line.strip().split("\t")
                records.append([sample, metric, status])

    df = pd.DataFrame(records, columns=["sample", "metric", "status"])
    return df

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fastq_dir", required=True)
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--threads", default=4, type=int)
    args = parser.parse_args()

    fastq_files = glob.glob(os.path.join(args.fastq_dir, "*.fastq.gz"))
    if len(fastq_files) == 0:
        raise FileNotFoundError("No FASTQ files found.")

    # 1) Run FastQC
    run_fastqc(fastq_files, args.out_dir, args.threads)

    # 2) Parse QC summaries
    qc_df = parse_fastqc_summary(args.out_dir)
    qc_df.to_csv(os.path.join(args.out_dir, "fastqc_summary.csv"), index=False)
    print(f"[INFO] QC summary saved: {args.out_dir}/fastqc_summary.csv")

if __name__ == "__main__":
    main()
