
# **README.md for `preprocessing/` **

# Preprocessing Pipeline

This directory contains the preprocessing components of the **ctDNA-analysis** workflow.  
It covers all steps from raw FASTQ files to variant-ready BAM and VCF outputs.  
The scripts are designed to be modular, reproducible, and suitable for both simulation-based and real ctDNA datasets.

---

## Overview of Steps

1. **FASTQ Quality Control**  
2. **BWA-MEM2 Alignment (with Read Group tags)**  
3. **UMI-Aware Deduplication (UMI-tools)**  
4. **Variant Calling (Mutect2 with optional PoN/Germline/Intervals)**  
5. **Pipeline execution: Linux shell or Snakemake**

Each component can be executed individually or as a unified pipeline.

---

## Files in This Directory

### **1. `01_quality_control.py`**  
Runs **FastQC** on raw FASTQ files and generates summary tables.  
- Input: `*.fastq.gz`  
- Output: FastQC reports + aggregated summary CSV  
- Purpose: basic QC before alignment.

---

### **2. `02_alignment.py`**  
Performs alignment using **BWA-MEM2**.  
Key features:
- Automatic detection of FASTQ pairs (`R1` / `R2`)  
- Read Group (`@RG`) tags added per sample  
- Sorted BAM generation  
- BAM indexing with `samtools index`  
- Assumes UMI extraction has already been performed.

Output example:  

results/alignment/sample.sorted.bam
results/alignment/sample.sorted.bam.bai


---

### **3. `03_umi_processing_umitools.py`**  
UMI-aware deduplication using **UMI-tools**.  
- Uses UMIs stored in the **RX tag**  
- Supports multiple deduplication methods (default: `directional`)  
- Multi-threaded UMI-tools + multi-threaded indexing  
- Generates metrics for each sample

Output example:  

results/umi_dedup/sample.dedup.bam
results/umi_dedup/sample.dedup.bam.bai
results/umi_dedup/sample.dedup_metrics.txt



---

### **4. `04_variant_calling_Mutect2.py`**  
Wrapper for **GATK Mutect2**, supporting both simulated and real ctDNA workflows.  
Features:
- Optional Panel of Normals (PoN)  
- Optional germline resource (e.g., gnomAD)  
- Optional interval BED file  
- Performs:
  1. Mutect2  
  2. LearnReadOrientationModel  
  3. FilterMutectCalls  

Output example:  

results/mutect/sample.vcf.gz
results/mutect/sample.unfiltered.vcf.gz
results/mutect/sample.f1r2.tar.gz
results/mutect/sample.orientation-model.tar.gz


---

### **5. `run_all.sh`**  
A full Linux-based pipeline script that sequentially runs:
1. QC  
2. Alignment  
3. UMI-tools dedup  
4. Mutect2 variant calling  

Highly useful for quick testing on HPC/Linux environments.

Usage:
```bash
bash run_all.sh
````

---

### **6. `Snakefile`**

Snakemake implementation of the entire preprocessing workflow.
Advantages:

* Automatic parallelization
* Built-in dependency management
* Easy reproducibility
* Suitable for production-scale pipelines

Run with:

```bash
snakemake -j 8
```

---

## Notes

* This preprocessing layer assumes that **UMI extraction** was completed upstream (e.g., `umi_tools extract`).
* All scripts follow a consistent directory structure under `../data/` and `../results/`.
* The preprocessing outputs feed directly into the downstream analysis scripts in the `analysis/` directory.

---

## Status

This preprocessing module is complete at Version 1 and will be extended with:

* MultiQC integration
* fgbio support (optional)
* Additional QC metrics for ultra-deep ctDNA


