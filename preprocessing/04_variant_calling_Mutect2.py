#!/usr/bin/env python3
"""
04_variant_calling.py
Mutect2 wrapper (GATK4) for ctDNA analysis.

Supports:
    - Optional Panel of Normals (PoN)
    - Optional germline resource (e.g., gnomAD)
    - Optional interval BED file
    - Ultra-deep ctDNA and simulated data

Inputs:
    - Deduplicated BAMs from 03_umi_processing_umitools.py
    - Reference genome (hg38 or hg19)

Usage example (simulation):
    python 04_variant_calling.py \
        --bam_dir ./results/umi_dedup \
        --reference ./refs/hg38.fa \
        --out_dir ./results/mutect

Usage example (real ctDNA):
    python 04_variant_calling.py \
        --bam_dir ./results/umi_dedup \
        --reference ./refs/hg38.fa \
        --out_dir ./results/mutect \
        --pon ./refs/pon.vcf.gz \
        --germline ./refs/af-only-gnomad.hg38.vcf.gz \
        --intervals ./refs/panel.bed \
        --threads 8
"""

import os
import argparse
import glob
import subprocess

def run_cmd(cmd):
    print("[CMD] " + " ".join(cmd))
    subprocess.run(cmd, check=True)

def find_bam_files(bam_dir):
    bams = sorted(glob.glob(os.path.join(bam_dir, "*.dedup.bam")))
    if len(bams) == 0:
        raise FileNotFoundError(f"No dedup.bam files found in: {bam_dir}")
    return bams

def run_mutect2(bams, reference, out_dir, pon, germline, intervals, threads):
    os.makedirs(out_dir, exist_ok=True)

    for bam in bams:
        sample = os.path.basename(bam).replace(".dedup.bam", "")
        vcf = os.path.join(out_dir, f"{sample}.vcf.gz")
        unfiltered_vcf = os.path.join(out_dir, f"{sample}.unfiltered.vcf.gz")
        f1r2 = os.path.join(out_dir, f"{sample}.f1r2.tar.gz")

        cmd = [
            "gatk", "Mutect2",
            "-R", reference,
            "-I", bam,
            "-O", unfiltered_vcf,
            "--f1r2-tar-gz", f1r2,
            "-nt", str(threads)  # legacy, GATK still accepts -nt/-nct
        ]

        if intervals:
            cmd += ["-L", intervals]

        if pon:
            cmd += ["--panel-of-normals", pon]

        if germline:
            cmd += ["--germline-resource", germline]

        print(f"[INFO] Running Mutect2 for sample: {sample}")
        run_cmd(cmd)

        # Learn read orientation model
        model = os.path.join(out_dir, f"{sample}.orientation-model.tar.gz")
        run_cmd([
            "gatk", "LearnReadOrientationModel",
            "-I", f1r2,
            "-O", model
        ])

        # Filter Mutect2 calls
        run_cmd([
            "gatk", "FilterMutectCalls",
            "-R", reference,
            "-V", unfiltered_vcf,
            "--orientation-bias-artifact-priors", model,
            "-O", vcf
        ])

        print(f"[INFO] Completed Mutect2: {vcf}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam_dir", required=True, help="Directory of .dedup.bam files")
    parser.add_argument("--reference", required=True, help="Reference FASTA")
    parser.add_argument("--out_dir", required=True, help="Output directory")
    parser.add_argument("--pon", default=None, help="Panel of Normals VCF (optional)")
    parser.add_argument("--germline", default=None, help="Germline resource (optional)")
    parser.add_argument("--intervals", default=None, help="BED or interval list (optional)")
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()

    bams = find_bam_files(args.bam_dir)

    run_mutect2(
        bams=bams,
        reference=args.reference,
        out_dir=args.out_dir,
        pon=args.pon,
        germline=args.germline,
        intervals=args.intervals,
        threads=args.threads
    )

if __name__ == "__main__":
    main()
