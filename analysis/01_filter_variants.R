#!/usr/bin/env Rscript

############################################################
# 01_filter_variants.R
# Purpose: Filter Mutect2 VCF results for ctDNA analysis
# Input:  Mutect2 VCF (sample.vcf.gz)
# Output: Filtered CSV table
#
# Notes:
# - ctDNA data often requires ultra-deep coverage
# - Filtering defaults here are conservative but adjustable
############################################################

suppressPackageStartupMessages({
  library(vcfR)
  library(dplyr)
  library(readr)
  library(optparse)
})

############################################################
# Command-line arguments
############################################################

option_list <- list(
  make_option(c("-i", "--vcf"), type="character", help="Input VCF file (.vcf.gz)", metavar="file"),
  make_option(c("-o", "--out"), type="character", help="Output CSV file", metavar="file"),
  make_option(c("--min_dp"), type="integer", default=1000,
              help="Minimum depth (default: 1000)"),
  make_option(c("--min_alt"), type="integer", default=3,
              help="Minimum alt count (default: 3)"),
  make_option(c("--min_vaf"), type="double", default=0.0005,
              help="Minimum VAF (default: 0.0005 = 0.05%)")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$vcf) | is.null(opt$out)) {
  stop("Missing required arguments: --vcf and --out")
}

############################################################
# Load VCF
############################################################

cat("Loading VCF:", opt$vcf, "\n")
vcf <- read.vcfR(opt$vcf)

fix <- vcf@fix
chrom  <- fix[, "CHROM"]
pos    <- as.integer(fix[, "POS"])
ref    <- fix[, "REF"]
alt    <- fix[, "ALT"]
filter <- fix[, "FILTER"]

############################################################
# Extract FORMAT fields
############################################################

dp  <- extract.gt(vcf, "DP", as.numeric = TRUE)
af  <- extract.gt(vcf, "AF", as.numeric = TRUE)
ad  <- extract.gt(vcf, "AD")

# Split AD into ref / alt counts
ad_split <- strsplit(ad, ",")
ref_count <- as.numeric(sapply(ad_split, `[`, 1))
alt_count <- as.numeric(sapply(ad_split, `[`, 2))

############################################################
# Combine into a data frame
############################################################

df <- data.frame(
  chrom, pos, ref, alt,
  dp, ref_count, alt_count, af, filter,
  stringsAsFactors = FALSE
)

############################################################
# Filtering rules (default ctDNA thresholds)
############################################################

filtered <- df %>%
  filter(filter == "PASS") %>%
  filter(!is.na(dp), !is.na(alt_count), !is.na(af)) %>%
  filter(dp >= opt$min_dp) %>%
  filter(alt_count >= opt$min_alt) %>%
  filter(af >= opt$min_vaf)

############################################################
# Save output
############################################################

cat("Writing filtered output to:", opt$out, "\n")
write_csv(filtered, opt$out)

cat("Done.\n")
