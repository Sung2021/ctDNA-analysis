#!/usr/bin/env Rscript

############################################################
# 03_visualization.R
# Purpose: Basic visualization of ctDNA variant profiles.
# Input:  Filtered variant table (CSV)
# Output: PNG figures (VAF distribution, depth-VAF, mutation spectrum)
#
# Notes:
# - Minimal but practical visualization set for ctDNA QC.
# - Designed for ultra-deep Mutect2 calls with VAF << 1%.
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(optparse)
})

############################################################
# Command-line arguments
############################################################

option_list <- list(
  make_option(c("-i", "--input"), type="character",
              help="Filtered variant CSV (output of 01_filter_variants.R)"),
  make_option(c("-o", "--outdir"), type="character", default="results/plots",
              help="Output directory for plots")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input)) {
  stop("Missing required argument: --input")
}

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

############################################################
# Load data
############################################################

cat("Loading:", opt$input, "\n")
df <- read_csv(opt$input, show_col_types = FALSE)

if (!all(c("chrom","pos","ref","alt","dp","alt_count","af") %in% colnames(df))) {
  stop("Input table missing required columns.")
}

############################################################
# 1) VAF distribution
############################################################

p1 <- ggplot(df, aes(x = af)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 0.01)) +
  labs(
    title = "VAF Distribution",
    x = "Variant Allele Frequency (VAF)",
    y = "Count"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(opt$outdir, "vaf_distribution.png"),
  plot = p1,
  width = 6, height = 4, dpi = 300
)

############################################################
# 2) Depth vs VAF
############################################################

p2 <- ggplot(df, aes(x = dp, y = af)) +
  geom_point(alpha = 0.6, color = "darkred") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.01)) +
  labs(
    title = "Depth vs VAF",
    x = "Depth (DP)",
    y = "VAF"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(opt$outdir, "depth_vs_vaf.png"),
  plot = p2,
  width = 6, height =
