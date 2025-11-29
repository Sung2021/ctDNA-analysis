# ctDNA-analysis

**Version 1 — 2025.11.29**

## Project Structure

```
ctDNA-analysis/
├── preprocessing/              
│   ├── 01_quality_control.py      # FASTQ QC
│   ├── 02_alignment.py            # BWA alignment
│   ├── 03_umi_processing.py       # UMI deduplication
│   └── 04_variant_calling.py      # MuTect2 wrapper
│
├── analysis/
│   ├── 01_filter_variants.R       # VAF filtering
│   ├── 02_annotation.R            # Variant annotation
│   ├── 03_visualization.R         # Result visualization
│   └── 04_statistics.R            # Statistical validation
│
├── simulation/
│   ├── generate_reads.py          # Synthetic read simulation
│   └── add_noise.R                # Noise modeling
│
├── data/                          # Raw / intermediate / processed data
└── docs/                          # Documentation, reports, references
```

