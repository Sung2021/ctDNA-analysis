
---

## Files in This Directory

### **1. `01_filter_variants.R`**
Filters Mutect2 VCF output into a high-confidence variant table.  
Main tasks:
- Extract depth (DP), allele depths (AD), VAF (AF), and FILTER status  
- Apply ctDNA thresholds (depth, alt count, VAF)  
- Output cleaned CSV ready for annotation

---

### **2. `02_annotation.R`**
Annotates variants using **VEP (Variant Effect Predictor)**.  
Features:
- Accepts VCF or filtered CSV  
- Adds functional consequence, gene symbol, transcript, SIFT/PolyPhen, and population AF  
- Outputs a transcript-level annotation table

---

### **3. `03_visualization.R`**
Generates essential QC plots:
- **VAF distribution**  
- **Depth vs VAF scatter**  
- **Mutation spectrum (REF>ALT)**  

These plots help validate variant calling and detect biases or anomalies.

---

### **4. `04_statistics.R`**
Evaluates variant-detection performance using **spike-in ground truth**.  
Provides:
- Detection sensitivity by expected VAF  
- Limit of Detection (LOD) curve  
- Per-variant detection status table  

Useful for benchmarking pipeline performance under ultra-deep ctDNA settings.

---

## Notes

- All scripts assume inputs from the preprocessing module.  
- Designed specifically for ultra-low VAF ctDNA conditions (<1%).  
- Spike-in simulation files (in `/simulation/`) integrate directly with `04_statistics.R`.  
- Outputs are stored in structured folders under `results/`.

---

## Status

The analysis module is functional for Version 1 and will be extended with:
- Multi-sample summarization  
- ROC/PR analysis for caller comparison  
- Additional visualization options  
