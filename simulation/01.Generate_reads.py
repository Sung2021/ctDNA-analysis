"""
ctDNA 변이콜링 시뮬레이션
대규모 데이터셋(6,000 포지션)을 생성하고 변이콜링 파이프라인을 실행합니다.
"""

import pandas as pd
import numpy as np
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

# ============================================================================
# 시뮬레이션 파라미터
# ============================================================================

# 주요 암 유전자의 실제 좌표 사용
target_regions = {
    'TP53': {'chr': 'chr17', 'start': 7571720, 'end': 7590868, 'positions_per_gene': 2500},
    'EGFR': {'chr': 'chr7', 'start': 55086714, 'end': 55275031, 'positions_per_gene': 2000},
    'KRAS': {'chr': 'chr12', 'start': 25357723, 'end': 25403870, 'positions_per_gene': 1500}
}

# UMI 및 PCR 파라미터
umi_mean = 500  # 평균 UMI 분자 수
umi_std = 100
avg_pcr_copies = 100  # PCR 평균 증폭 횟수
pcr_std = 30

# VAF 분포 정의
tumor_fractions = {
    'high': 0.05,      # 5% - 명확한 변이
    'medium': 0.015,   # 1.5% - 감지 가능
    'low': 0.005       # 0.5% - 경계 케이스
}

# 에러율
sequencing_error_rate = 0.0005  # Sequencing error rate
umi_error_rate = 0.0001  # UMI error rate

# Hotspot 변이 정의 (유전자별)
hotspot_variants = {
    'TP53': [
        {'pos': 7577548, 'vaf': 'high', 'variant': 'R248Q'},
        {'pos': 7577559, 'vaf': 'high', 'variant': 'R248W'},
        {'pos': 7578190, 'vaf': 'medium', 'variant': 'R175H'},
        {'pos': 7578191, 'vaf': 'low', 'variant': 'Y234H'},
    ],
    'EGFR': [
        {'pos': 55191822, 'vaf': 'high', 'variant': 'L858R'},
        {'pos': 55191817, 'vaf': 'medium', 'variant': 'L861Q'},
        {'pos': 55086776, 'vaf': 'medium', 'variant': 'E746-A750'},
    ],
    'KRAS': [
        {'pos': 25398284, 'vaf': 'high', 'variant': 'G12C'},
        {'pos': 25398285, 'vaf': 'high', 'variant': 'G12V'},
        {'pos': 25398286, 'vaf': 'medium', 'variant': 'G13D'},
        {'pos': 25380275, 'vaf': 'low', 'variant': 'Q61H'},
    ]
}

# ============================================================================
# 함수 정의
# ============================================================================

def generate_umi_count(mean=umi_mean, std=umi_std):
    """UMI 분자 수 생성 (정규분포)"""
    return max(int(np.random.normal(mean, std)), 10)


def binomial_test_pvalue(alt_reads, total_reads, error_rate):
    """
    이항분포 검정으로 p-value 계산
    NGS 변이콜링의 기본 통계적 검증 방법
    """
    if total_reads < 2:
        return 1.0
    return 1 - stats.binom.cdf(alt_reads - 1, total_reads, error_rate)


def _generate_molecule_counts(n_unique_molecules, true_vaf):
    """UMI 분자 수와 true_vaf를 기반으로 변이/정상 분자 수를 결정"""
    n_variant_molecules = np.random.binomial(n_unique_molecules, true_vaf)
    n_ref_molecules = n_unique_molecules - n_variant_molecules
    return n_variant_molecules, n_ref_molecules


def _simulate_pcr_amplification(n_molecules):
    """주어진 분자 수에 대해 PCR 증폭 후 총 리드 수를 반환"""
    total_reads = 0
    for _ in range(n_molecules):
        copies = max(1, int(np.random.normal(avg_pcr_copies, pcr_std)))
        total_reads += copies
    return total_reads


def generate_variant_data(pos, gene_name, region, hotspots):
    """
    단일 유전체 위치에 대한 시뮬레이션 데이터를 생성
    (통계 분석 및 필터링은 포함하지 않음)
    """
    n_unique_molecules = generate_umi_count()
    
    # 변이 정보 초기화
    n_variant_molecules, n_ref_molecules = 0, n_unique_molecules
    variant_reads, ref_reads, total_depth = 0, 0, 0
    is_true_variant = False
    variant_type = 'background'
    true_vaf = None

    if pos in hotspots:
        # ===== 🎯 진짜 변이 (Hotspot) 시뮬레이션 =====
        hotspot_info = hotspots[pos]
        true_vaf = tumor_fractions[hotspot_info['vaf']]
        
        # UMI 수준 변이 결정
        n_variant_molecules, n_ref_molecules = _generate_molecule_counts(n_unique_molecules, true_vaf)
        
        # PCR 증폭 시뮬레이션
        variant_reads = _simulate_pcr_amplification(n_variant_molecules)
        ref_reads = _simulate_pcr_amplification(n_ref_molecules)
        total_depth = variant_reads + ref_reads
        
        is_true_variant = True
        variant_type = hotspot_info['variant']
        
    else:
        # ===== 👻 배경 잡음 시뮬레이션 (Sequencing Error) =====
        total_depth = _simulate_pcr_amplification(n_unique_molecules)
        variant_reads = np.random.binomial(total_depth, sequencing_error_rate)
        ref_reads = total_depth - variant_reads
    
    vaf = variant_reads / total_depth if total_depth > 0 else 0
    
    return {
        'chrom': region['chr'],
        'position': pos,
        'gene': gene_name,
        'variant_type': variant_type,
        'is_true_variant': is_true_variant,
        'true_vaf': true_vaf,
        
        'n_unique_molecules': n_unique_molecules,
        'n_alt_molecules': n_variant_molecules,
        'n_ref_molecules': n_ref_molecules,
        
        'alt_reads': variant_reads,
        'ref_reads': ref_reads,
        'total_reads': total_depth,
        'vaf': vaf,
        'vaf_percent': vaf * 100
    }


def analyze_and_filter_variant(data_record):
    """생성된 데이터 레코드에 통계 검정 및 필터 기준을 적용"""
    
    variant_reads = data_record['alt_reads']
    total_depth = data_record['total_reads']
    vaf = data_record['vaf']
    n_variant_molecules = data_record['n_alt_molecules']
    
    # 통계 검정
    binomial_pval = binomial_test_pvalue(
        variant_reads, 
        total_reads=total_depth, 
        error_rate=sequencing_error_rate
    )
    
    log10_pval = np.log10(binomial_pval) if binomial_pval > 0 else -300
    
    # ===== 📏 필터링 로직 =====
    
    # Stringent (엄격한) 필터 기준
    pass_stringent = (
        variant_reads >= 5 and
        n_variant_molecules >= 3 and
        vaf >= 0.005 and
        binomial_pval < 0.01 and
        total_depth >= 10000
    )
    
    # Sensitive (민감한) 필터 기준
    pass_sensitive = (
        variant_reads >= 3 and
        n_variant_molecules >= 2 and
        vaf >= 0.001 and
        binomial_pval < 0.05 and
        total_depth >= 5000
    )
    
    # 결과 레코드 업데이트
    data_record.update({
        'binomial_pval': binomial_pval,
        'log10_pval': log10_pval,
        'pass_stringent': pass_stringent,
        'pass_sensitive': pass_sensitive,
        
        # 개별 필터 상태
        'pass_depth_stringent': total_depth >= 10000,
        'pass_depth_sensitive': total_depth >= 5000,
        'pass_alt_reads_stringent': variant_reads >= 5,
        'pass_alt_reads_sensitive': variant_reads >= 3,
        'pass_vaf_stringent': vaf >= 0.005,
        'pass_vaf_sensitive': vaf >= 0.001,
        'pass_pval_stringent': binomial_pval < 0.01,
        'pass_pval_sensitive': binomial_pval < 0.05,
        'pass_umi_stringent': n_variant_molecules >= 3,
        'pass_umi_sensitive': n_variant_molecules >= 2,
    })
    
    return data_record


def simulate_gene_variants(gene_name, region, n_positions):
    """
    특정 유전자의 변이 시뮬레이션을 실행하고, 각 포지션에 대해 분석 및 필터링을 적용
    """
    hotspots_dict = {v['pos']: v for v in hotspot_variants.get(gene_name, [])}
    hotspot_positions = set(hotspots_dict.keys())
    
    # 배경 포지션 생성 (핫스팟 위치 제외)
    all_uniform_positions = set(np.linspace(region['start'], region['end'], n_positions, dtype=int))
    background_positions = list(all_uniform_positions - hotspot_positions)
    
    # 핫스팟 위치와 배경 위치 통합
    simulation_positions = list(hotspot_positions) + background_positions
    
    data = []
    
    for pos in simulation_positions:
        data_record = generate_variant_data(pos, gene_name, region, hotspots_dict)
        final_record = analyze_and_filter_variant(data_record)
        data.append(final_record)
    
    return data


# ============================================================================
# 메인 시뮬레이션 실행
# ============================================================================

if __name__ == "__main__":
    print("=" * 80)
    print("ctDNA 변이콜링 시뮬레이션 시작")
    print("=" * 80)

    all_data = []

    for gene_name, region in target_regions.items():
        print(f"\n[{gene_name}] 시뮬레이션 중...")
        gene_data = simulate_gene_variants(
            gene_name, 
            region, 
            region['positions_per_gene']
        )
        all_data.extend(gene_data)
        print(f"  - {len(gene_data)} 포지션 생성 완료")

    sim_ctdna = pd.DataFrame(all_data)

    # ============================================================================
    # 결과 분석 및 출력
    # ============================================================================

    print("\n" + "=" * 80)
    print("시뮬레이션 요약")
    print("=" * 80)

    print(f"\n📊 데이터셋 규모:")
    print(f"  - 전체 포지션: {len(sim_ctdna):,}")
    print(f"  - TP53: {len(sim_ctdna[sim_ctdna['gene']=='TP53']):,}")
    print(f"  - EGFR: {len(sim_ctdna[sim_ctdna['gene']=='EGFR']):,}")
    print(f"  - KRAS: {len(sim_ctdna[sim_ctdna['gene']=='KRAS']):,}")

    print(f"\n📈 시퀀싱 통계:")
    print(f"  - 평균 깊이: {sim_ctdna['total_reads'].mean():,.0f}")
    print(f"  - 깊이 범위: {sim_ctdna['total_reads'].min():,} - {sim_ctdna['total_reads'].max():,}")
    print(f"  - 평균 UMI 분자: {sim_ctdna['n_unique_molecules'].mean():.1f}")

    print(f"\n🔍 변이 정보:")
    print(f"  - 진정 변이 (True variant): {sim_ctdna['is_true_variant'].sum()}")
    print(f"  - 배경 노이즈: {(~sim_ctdna['is_true_variant']).sum()}")

    # 두 가지 필터 결과 비교
    print(f"\n✅ 필터링 결과 비교:")
    print(f"\n  [Stringent Filter]")
    stringent_total = sim_ctdna['pass_stringent'].sum()
    print(f"  - 통과: {stringent_total} ({stringent_total/len(sim_ctdna)*100:.2f}%)")

    print(f"\n  [Sensitive Filter]")
    sensitive_total = sim_ctdna['pass_sensitive'].sum()
    print(f"  - 통과: {sensitive_total} ({sensitive_total/len(sim_ctdna)*100:.2f}%)")

    print(f"\n📋 진정 변이 상세:")
    true_vars = sim_ctdna[sim_ctdna['is_true_variant']].sort_values('vaf_percent', ascending=False)
    print(true_vars[['gene', 'position', 'variant_type', 'n_alt_molecules', 'alt_reads', 
                     'total_reads', 'vaf_percent', 'pass_stringent', 'pass_sensitive']].to_string(index=False))

    # ============================================================================
    # 📊 성능 지표 계산 및 출력
    # ============================================================================

    true_variants = sim_ctdna[sim_ctdna['is_true_variant'] == True]
    background_noise = sim_ctdna[sim_ctdna['is_true_variant'] == False]

    print("\n" + "=" * 80)
    print("✅ NGS 시스템 성능 평가")
    print("=" * 80)

    for filter_name in ['stringent', 'sensitive']:
        print(f"\n--- [{filter_name.capitalize()} Filter] ---")
        
        tp_count = true_variants[f'pass_{filter_name}'].sum()
        sensitivity = tp_count / len(true_variants) if len(true_variants) > 0 else 0
        
        print(f"  - 🎯 민감도 (Sensitivity): {sensitivity * 100:.2f}% ({tp_count}/{len(true_variants)} 변이 검출)")
        
        fp_count = background_noise[f'pass_{filter_name}'].sum()
        fpr = fp_count / len(background_noise) if len(background_noise) > 0 else 0
        
        print(f"  - 👻 위양성률 (FPR): {fpr * 100:.4f}% ({fp_count}/{len(background_noise)} 노이즈 오인)")
        
        low_vaf_success = true_variants[
            (true_variants['true_vaf'] == 0.005) & (true_variants[f'pass_{filter_name}'])
        ]
        print(f"  - Low VAF (0.5%) 검출: {len(low_vaf_success)}/{len(true_variants[true_variants['true_vaf'] == 0.005])}개")

    # ============================================================================
    # 파일 저장
    # ============================================================================

    output_file = 'ctdna_simulated_data.csv'
    sim_ctdna.to_csv(output_file, index=False)
    print(f"\n💾 데이터 저장: {output_file}")

    called_variants_stringent = sim_ctdna[sim_ctdna['pass_stringent']].copy()
    called_variants_stringent_file = 'ctdna_called_variants_stringent.csv'
    called_variants_stringent.to_csv(called_variants_stringent_file, index=False)
    print(f"📌 Stringent 변이콜링 결과: {called_variants_stringent_file}")

    called_variants_sensitive = sim_ctdna[sim_ctdna['pass_sensitive']].copy()
    called_variants_sensitive_file = 'ctdna_called_variants_sensitive.csv'
    called_variants_sensitive.to_csv(called_variants_sensitive_file, index=False)
    print(f"📌 Sensitive 변이콜링 결과: {called_variants_sensitive_file}")

    vaf_summary = pd.DataFrame({
        'gene': sim_ctdna['gene'],
        'vaf': sim_ctdna['vaf'],
        'vaf_percent': sim_ctdna['vaf_percent'],
        'is_true': sim_ctdna['is_true_variant'],
        'pass_stringent': sim_ctdna['pass_stringent'],
        'pass_sensitive': sim_ctdna['pass_sensitive']
    })
    vaf_summary.to_csv('ctdna_vaf_distribution.csv', index=False)
    print(f"📊 VAF 분포: ctdna_vaf_distribution.csv")
    
    print("\n" + "=" * 80)
    print("시뮬레이션 완료!")
    print("=" * 80)
