# ctDNA-analysis

(Version1 on 25.11.29)

ctDNA-analysis/
├── preprocessing/
│   ├── 01_quality_control.py      # FASTQ QC
│   ├── 02_alignment.py            # BWA 실행
│   ├── 03_umi_processing.py       # UMI deduplication
│   └── 04_variant_calling.py      # MuTect2 wrapper
├── analysis/
│   ├── 01_filter_variants.R       # VAF 필터링
│   ├── 02_annotation.R            # 변이 주석
│   ├── 03_visualization.R         # 결과 시각화
│   └── 04_statistics.R            # 통계 검증
├── simulation/
│   ├── generate_reads.py          # 시뮬레이션 데이터
│   └── add_noise.R                # 노이즈 모델링
├── data/
└── docs/
