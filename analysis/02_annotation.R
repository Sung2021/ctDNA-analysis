#!/usr/bin/env Rscript

############################################################
# 02_annotation.R
# Purpose: Annotate variants using VEP output
# Input:  VCF (Mutect2 output) or filtered CSV
# Output: Annotated TSV using VEP annotation fields
#
# Notes:
# - Designed for ctDNA variants (low VAF)
# - Requires VEP pre-run using CLI: vep --vcf --everything
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(optparse)
})

############################################################
# Command-line arguments
############################################################

option_list <- list(
  make_option(c("-i", "--input"), type="character",
              help="Input file (VCF or CSV)"),
  make_option(c("-t", "--type"), type="character", default="vcf",
              help="Input type: vcf or csv (default: vcf)"),
  make_option(c("-a", "--annot"), type="character",
              help="VEP output file (annotated VCF or TXT)"),
  make_option(c("-o", "--out"), type="character",
              help="Output TSV file")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input) | is.null(opt$annot) | is.null(opt$out)) {
  stop("Missing required arguments: --input, --annot, --out")
}

############################################################
# Load input variant list
############################################################

cat("Loading input:", opt$input, "\n")

if (opt$type == "csv") {
  variants <- read_csv(opt$input, show_col_types = FALSE)
} else if (opt$type == "vcf") {
  library(vcfR)
  v <- read.vcfR(opt$input)

  variants <- data.frame(
    chrom = v@fix[, "CHROM"],
    pos   = as.integer(v@fix[, "POS"]),
    ref   = v@fix[, "REF"],
    alt   = v@fix[, "ALT"],
    stringsAsFactors = FALSE
  )
} else {
  stop("Unknown input type: use 'csv' or 'vcf'")
}

############################################################
# Load VEP annotations
############################################################

cat("Loading VEP annotation:", opt$annot, "\n")

# VEP standard output: one line per transcript consequence
annot_raw <- read_tsv(opt$annot, comment = "#", show_col_types = FALSE)

# VEP key fields
wanted_cols <- c(
  "Location", "Allele", "Consequence", "IMPACT",
  "SYMBOL", "Gene", "Feature_type", "Feature",
  "HGVSc", "HGVSp", "EXON", "INTRON",
  "Existing_variation", "AF", "gnomAD_AF",
  "SIFT", "PolyPhen", "CLIN_SIG"
)

annot <- annot_raw %>%
  select(any_of(wanted_cols)) %>%
  distinct()

############################################################
# Merge input variants with annotation
############################################################

cat("Merging variants with annotation...\n")

# VEP "Location" field looks like: "1:12345"
annot_split <- annot %>%
  mutate(
    chrom = sub(":.*", "", Location),
    pos   = as.integer(sub(".*:", "", Location))
  )

out <- variants %>%
  left_join(annot_split, by = c("chrom", "pos"))

############################################################
# Save result
############################################################

cat("Writing output:", opt$out, "\n")
write_tsv(out, opt$out)

cat("Done.\n")
