# ============================================================================
# ctDNA 변이콜링 시뮬레이션 (R 버전)
# ============================================================================

library(dplyr)
library(tidyr)

set.seed(42)

# ============================================================================
# 시뮬레이션 파라미터
# ============================================================================

# 주요 암 유전자의 실제 좌표
target_regions <- list(
  TP53 = list(chr = "chr17", start = 7571720, end = 7590868, n_positions = 2500),
  EGFR = list(chr = "chr7", start = 55086714, end = 55275031, n_positions = 2000),
  KRAS = list(chr = "chr12", start = 25357723, end = 25403870, n_positions = 1500)
)

# UMI 및 PCR 파라미터
umi_mean <- 500
umi_std <- 100
avg_pcr_copies <- 100
pcr_std <- 30

# VAF 분포 정의
tumor_fractions <- list(
  high = 0.05,
  medium = 0.015,
  low = 0.005
)

# 에러율
sequencing_error_rate <- 0.0005
umi_error_rate <- 0.0001

# Hotspot 변이 정의
hotspot_variants <- list(
  TP53 = data.frame(
    pos = c(7577548, 7577559, 7578190, 7578191),
    vaf = c("high", "high", "medium", "low"),
    variant = c("R248Q", "R248W", "R175H", "Y234H"),
    stringsAsFactors = FALSE
  ),
  EGFR = data.frame(
    pos = c(55191822, 55191817, 55086776),
    vaf = c("high", "medium", "medium"),
    variant = c("L858R", "L861Q", "E746-A750"),
    stringsAsFactors = FALSE
  ),
  KRAS = data.frame(
    pos = c(25398284, 25398285, 25398286, 25380275),
    vaf = c("high", "high", "medium", "low"),
    variant = c("G12C", "G12V", "G13D", "Q61H"),
    stringsAsFactors = FALSE
  )
)

# ============================================================================
# 함수 정의
# ============================================================================

generate_umi_count <- function(mean = umi_mean, std = umi_std) {
  # UMI 분자 수 생성
  max(as.integer(rnorm(1, mean, std)), 10)
}

binomial_test_pvalue <- function(alt_reads, total_reads, error_rate) {
  # 이항분포 검정
  if (total_reads < 2) return(1.0)
  1 - pbinom(alt_reads - 1, total_reads, error_rate)
}

generate_molecule_counts <- function(n_unique_molecules, true_vaf) {
  # UMI 분자 수와 true_vaf 기반으로 변이/정상 분자 수 결정
  n_variant_molecules <- rbinom(1, n_unique_molecules, true_vaf)
  n_ref_molecules <- n_unique_molecules - n_variant_molecules
  list(variant = n_variant_molecules, ref = n_ref_molecules)
}

simulate_pcr_amplification <- function(n_molecules) {
  # PCR 증폭 시뮬레이션
  if (n_molecules == 0) return(0)
  
  total_reads <- 0
  for (i in 1:n_molecules) {
    copies <- max(1, as.integer(rnorm(1, avg_pcr_copies, pcr_std)))
    total_reads <- total_reads + copies
  }
  total_reads
}

generate_variant_data <- function(pos, gene_name, region, hotspots_df) {
  # 단일 유전체 위치에 대한 시뮬레이션 데이터 생성
  
  n_unique_molecules <- generate_umi_count()
  
  # 초기값
  n_variant_molecules <- 0
  n_ref_molecules <- n_unique_molecules
  variant_reads <- 0
  ref_reads <- 0
  total_depth <- 0
  is_true_variant <- FALSE
  variant_type <- "background"
  true_vaf <- NA
  
  # Hotspot 확인
  is_hotspot <- pos %in% hotspots_df$pos
  
  if (is_hotspot) {
    # 진짜 변이
    hotspot_info <- hotspots_df[hotspots_df$pos == pos, ]
    vaf_category <- hotspot_info$vaf[1]
    true_vaf <- tumor_fractions[[vaf_category]]
    
    # UMI 수준 변이 결정
    molecule_counts <- generate_molecule_counts(n_unique_molecules, true_vaf)
    n_variant_molecules <- molecule_counts$variant
    n_ref_molecules <- molecule_counts$ref
    
    # PCR 증폭
    variant_reads <- simulate_pcr_amplification(n_variant_molecules)
    ref_reads <- simulate_pcr_amplification(n_ref_molecules)
    total_depth <- variant_reads + ref_reads
    
    is_true_variant <- TRUE
    variant_type <- hotspot_info$variant[1]
    
  } else {
    # 배경 노이즈
    total_depth <- simulate_pcr_amplification(n_unique_molecules)
    variant_reads <- rbinom(1, total_depth, sequencing_error_rate)
    ref_reads <- total_depth - variant_reads
  }
  
  vaf <- ifelse(total_depth > 0, variant_reads / total_depth, 0)
  
  data.frame(
    chrom = region$chr,
    position = pos,
    gene = gene_name,
    variant_type = variant_type,
    is_true_variant = is_true_variant,
    true_vaf = true_vaf,
    n_unique_molecules = n_unique_molecules,
    n_alt_molecules = n_variant_molecules,
    n_ref_molecules = n_ref_molecules,
    alt_reads = variant_reads,
    ref_reads = ref_reads,
    total_reads = total_depth,
    vaf = vaf,
    vaf_percent = vaf * 100,
    stringsAsFactors = FALSE
  )
}

analyze_and_filter_variant <- function(data_record) {
  # 통계 검정 및 필터 기준 적용
  
  variant_reads <- data_record$alt_reads
  total_depth <- data_record$total_reads
  vaf <- data_record$vaf
  n_variant_molecules <- data_record$n_alt_molecules
  
  # 통계 검정
  binomial_pval <- binomial_test_pvalue(
    variant_reads, 
    total_depth, 
    sequencing_error_rate
  )
  
  log10_pval <- ifelse(binomial_pval > 0, log10(binomial_pval), -300)
  
  # Stringent 필터
  pass_stringent <- (
    variant_reads >= 5 &
    n_variant_molecules >= 3 &
    vaf >= 0.005 &
    binomial_pval < 0.01 &
    total_depth >= 10000
  )
  
  # Sensitive 필터
  pass_sensitive <- (
    variant_reads >= 3 &
    n_variant_molecules >= 2 &
    vaf >= 0.001 &
    binomial_pval < 0.05 &
    total_depth >= 5000
  )
  
  # 결과에 추가
  data_record$binomial_pval <- binomial_pval
  data_record$log10_pval <- log10_pval
  data_record$pass_stringent <- pass_stringent
  data_record$pass_sensitive <- pass_sensitive
  
  # 개별 필터 상태
  data_record$pass_depth_stringent <- total_depth >= 10000
  data_record$pass_depth_sensitive <- total_depth >= 5000
  data_record$pass_alt_reads_stringent <- variant_reads >= 5
  data_record$pass_alt_reads_sensitive <- variant_reads >= 3
  data_record$pass_vaf_stringent <- vaf >= 0.005
  data_record$pass_vaf_sensitive <- vaf >= 0.001
  data_record$pass_pval_stringent <- binomial_pval < 0.01
  data_record$pass_pval_sensitive <- binomial_pval < 0.05
  data_record$pass_umi_stringent <- n_variant_molecules >= 3
  data_record$pass_umi_sensitive <- n_variant_molecules >= 2
  
  data_record
}

simulate_gene_variants <- function(gene_name, region, n_positions) {
  # 특정 유전자의 변이 시뮬레이션
  
  cat(sprintf("\n[%s] 시뮬레이션 중...\n", gene_name))
  
  hotspots_df <- hotspot_variants[[gene_name]]
  hotspot_positions <- hotspots_df$pos
  
  # 배경 포지션 생성
  all_positions <- as.integer(seq(region$start, region$end, length.out = n_positions))
  background_positions <- setdiff(all_positions, hotspot_positions)
  
  # 통합
  simulation_positions <- c(hotspot_positions, background_positions)
  
  # 각 포지션 시뮬레이션
  results_list <- lapply(simulation_positions, function(pos) {
    data_record <- generate_variant_data(pos, gene_name, region, hotspots_df)
    analyze_and_filter_variant(data_record)
  })
  
  # 데이터프레임으로 결합
  do.call(rbind, results_list)
}

# ============================================================================
# 메인 시뮬레이션 실행
# ============================================================================

cat(rep("=", 80), "\n", sep = "")
cat("ctDNA 변이콜링 시뮬레이션 시작 (R 버전)\n")
cat(rep("=", 80), "\n", sep = "")

all_data_list <- list()

for (gene_name in names(target_regions)) {
  region <- target_regions[[gene_name]]
  gene_data <- simulate_gene_variants(
    gene_name, 
    region, 
    region$n_positions
  )
  all_data_list[[gene_name]] <- gene_data
  cat(sprintf("  - %d 포지션 생성 완료\n", nrow(gene_data)))
}

sim_ctdna <- do.call(rbind, all_data_list)
rownames(sim_ctdna) <- NULL

# ============================================================================
# 결과 분석 및 출력
# ============================================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("시뮬레이션 요약\n")
cat(rep("=", 80), "\n", sep = "")

cat(sprintf("\n📊 데이터셋 규모:\n"))
cat(sprintf("  - 전체 포지션: %s\n", format(nrow(sim_ctdna), big.mark = ",")))
cat(sprintf("  - TP53: %s\n", format(sum(sim_ctdna$gene == "TP53"), big.mark = ",")))
cat(sprintf("  - EGFR: %s\n", format(sum(sim_ctdna$gene == "EGFR"), big.mark = ",")))
cat(sprintf("  - KRAS: %s\n", format(sum(sim_ctdna$gene == "KRAS"), big.mark = ",")))

cat(sprintf("\n📈 시퀀싱 통계:\n"))
cat(sprintf("  - 평균 깊이: %s\n", format(round(mean(sim_ctdna$total_reads)), big.mark = ",")))
cat(sprintf("  - 깊이 범위: %s - %s\n", 
            format(min(sim_ctdna$total_reads), big.mark = ","),
            format(max(sim_ctdna$total_reads), big.mark = ",")))
cat(sprintf("  - 평균 UMI 분자: %.1f\n", mean(sim_ctdna$n_unique_molecules)))

cat(sprintf("\n🔍 변이 정보:\n"))
cat(sprintf("  - 진정 변이 (True variant): %d\n", sum(sim_ctdna$is_true_variant)))
cat(sprintf("  - 배경 노이즈: %d\n", sum(!sim_ctdna$is_true_variant)))

cat(sprintf("\n✅ 필터링 결과 비교:\n"))
cat(sprintf("\n  [Stringent Filter]\n"))
stringent_total <- sum(sim_ctdna$pass_stringent)
cat(sprintf("  - 통과: %d (%.2f%%)\n", stringent_total, 
            stringent_total / nrow(sim_ctdna) * 100))

cat(sprintf("\n  [Sensitive Filter]\n"))
sensitive_total <- sum(sim_ctdna$pass_sensitive)
cat(sprintf("  - 통과: %d (%.2f%%)\n", sensitive_total,
            sensitive_total / nrow(sim_ctdna) * 100))

cat(sprintf("\n📋 진정 변이 상세:\n"))
true_vars <- sim_ctdna %>%
  filter(is_true_variant == TRUE) %>%
  arrange(desc(vaf_percent)) %>%
  select(gene, position, variant_type, n_alt_molecules, alt_reads, 
         total_reads, vaf_percent, pass_stringent, pass_sensitive)

print(true_vars, row.names = FALSE)

# ============================================================================
# 성능 지표 계산
# ============================================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("✅ NGS 시스템 성능 평가\n")
cat(rep("=", 80), "\n", sep = "")

true_variants <- sim_ctdna %>% filter(is_true_variant == TRUE)
background_noise <- sim_ctdna %>% filter(is_true_variant == FALSE)

for (filter_name in c("stringent", "sensitive")) {
  cat(sprintf("\n--- [%s Filter] ---\n", 
              tools::toTitleCase(filter_name)))
  
  filter_col <- paste0("pass_", filter_name)
  
  # 민감도
  tp_count <- sum(true_variants[[filter_col]])
  sensitivity <- tp_count / nrow(true_variants)
  cat(sprintf("  - 🎯 민감도 (Sensitivity): %.2f%% (%d/%d 변이 검출)\n",
              sensitivity * 100, tp_count, nrow(true_variants)))
  
  # 위양성률
  fp_count <- sum(background_noise[[filter_col]])
  fpr <- fp_count / nrow(background_noise)
  cat(sprintf("  - 👻 위양성률 (FPR): %.4f%% (%d/%d 노이즈 오인)\n",
              fpr * 100, fp_count, nrow(background_noise)))
  
  # Low VAF 검출
  low_vaf_variants <- true_variants %>% filter(true_vaf == 0.005)
  low_vaf_success <- sum(low_vaf_variants[[filter_col]])
  cat(sprintf("  - Low VAF (0.5%%) 검출: %d/%d개\n",
              low_vaf_success, nrow(low_vaf_variants)))
}

# ============================================================================
# 파일 저장
# ============================================================================

write.csv(sim_ctdna, "ctdna_simulated_data.csv", row.names = FALSE)
cat(sprintf("\n💾 데이터 저장: ctdna_simulated_data.csv\n"))

stringent_variants <- sim_ctdna %>% filter(pass_stringent == TRUE)
write.csv(stringent_variants, "ctdna_called_variants_stringent.csv", row.names = FALSE)
cat(sprintf("📌 Stringent 변이콜링 결과: ctdna_called_variants_stringent.csv\n"))

sensitive_variants <- sim_ctdna %>% filter(pass_sensitive == TRUE)
write.csv(sensitive_variants, "ctdna_called_variants_sensitive.csv", row.names = FALSE)
cat(sprintf("📌 Sensitive 변이콜링 결과: ctdna_called_variants_sensitive.csv\n"))

vaf_summary <- sim_ctdna %>%
  select(gene, vaf, vaf_percent, is_true_variant, pass_stringent, pass_sensitive) %>%
  rename(is_true = is_true_variant)

write.csv(vaf_summary, "ctdna_vaf_distribution.csv", row.names = FALSE)
cat(sprintf("📊 VAF 분포: ctdna_vaf_distribution.csv\n"))

cat("\n", rep("=", 80), "\n", sep = "")
cat("시뮬레이션 완료!\n")
cat(rep("=", 80), "\n", sep = "")