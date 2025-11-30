#!/usr/bin/env Rscript

############################################################
# 04_statistics.R
# Purpose:
#   Evaluate ctDNA variant detection performance using
#   spike-in simulation ground truth.
#
# Inputs:
#   - filtered variant CSV  (01_filter_variants.R output)
#   - spike-in truth table   (expected variants + expected VAF)
#
# Outputs:
#   - sensitivity table
#   - LOD curve plot
#   - detection summary
#
# Notes:
#   This script assumes ultra-deep ctDNA data where VAF < 1%.
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(optparse)
})

############################################################
# Command-line arguments
############################################################

option_list <- list(
  make_option(c("-d", "--detected"), type="character",
              help="Filtered variant CSV (detected variants)"),
  make_option(c("-t", "--truth"), type="character",
              help="Spike-in truth table (chrom,pos,ref,alt,expected_vaf)"),
  make_option(c("-o", "--outdir"), type="character", default="results/statistics",
              help="Output directory")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$detected) | is.null(opt$truth)) {
  stop("Required arguments: --detected and --truth")
}

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

############################################################
# Load data
############################################################

cat("Loading detected variants:", opt$detected, "\n")
detected <- read_csv(opt$detected, show_col_types = FALSE)

cat("Loading truth table:", opt$truth, "\n")
truth <- read_csv(opt$truth, show_col_types = FALSE)

# Expected columns:
# truth: chrom, pos, ref, alt, expected_vaf

############################################################
# Standardize types
############################################################

detected <- detected %>%
  mutate(pos = as.integer(pos))

truth <- truth %>%
  mutate(pos = as.integer(pos))

############################################################
# Match detected variants with truth (True Positives)
############################################################

matched <- truth %>%
  left_join(
    detected %>% select(chrom, pos, ref, alt, af),
    by = c("chrom", "pos", "ref", "alt")
  ) %>%
  mutate(
    detected_flag = ifelse(is.na(af), 0, 1)
  )

############################################################
# Sensitivity by VAF bin
############################################################

# Define bins for expected VAF
matched <- matched %>%
  mutate(
    vaf_bin = cut(
      expected_vaf,
      breaks = c(0, 0.0005, 0.001, 0.005, 0.01, 0.05, 1),
      labels = c("<0.05%", "0.05–0.1%", "0.1–0.5%", "0.5–1%", "1–5%", ">5%"),
      include.lowest = TRUE
    )
  )

sensitivity_tbl <- matched %>%
  group_by(vaf_bin) %>%
  summarise(
    total = n(),
    detected = sum(detected_flag),
    sensitivity = detected / total
  ) %>%
  ungroup()

write_csv(sensitivity_tbl, file.path(opt$outdir, "sensitivity_by_vaf.csv"))

############################################################
# LOD curve (expected VAF vs. detection probability)
############################################################

p1 <- ggplot(sensitivity_tbl, aes(x = vaf_bin, y = sensitivity)) +
  geom_col(fill = "steelblue") +
  ylim(0, 1) +
  labs(
    title = "Sensitivity vs Expected VAF (LOD Curve)",
    x = "Expected VAF (binned)",
    y = "Detection Rate"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(opt$outdir, "lod_curve.png"),
  plot = p1,
  width = 6, height = 4, dpi = 300
)

############################################################
# Save matched table (per-variant detection status)
############################################################

write_csv(matched, file.path(opt$outdir, "variant_detection_status.csv"))

cat("Completed. Results saved in:", opt$outdir, "\n")
