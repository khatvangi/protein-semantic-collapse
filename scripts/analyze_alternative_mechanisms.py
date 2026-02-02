#!/usr/bin/env python3
"""
Analyze alternative regulatory mechanisms in RBD families that DON'T use D/E autoinhibition.

Compares:
- D/E-enriched families: RRM_1, LSM, DEAD
- Non-D/E families: KH_1, S1, cold-shock (CSD), dsRBD

Looks for:
1. RGG/RG motif enrichment (arginine-glycine repeats)
2. SR-rich regions (serine-arginine, splicing related)
3. PQ-rich regions (proline-glutamine, aggregation prone)
4. Different IDR amino acid compositions
5. PTM site density differences
"""

import json
import re
from collections import defaultdict, Counter
from pathlib import Path
import numpy as np
from scipy import stats

# output directory
RESULTS_DIR = Path("/storage/kiran-stuff/protein-semantic-collapse/results")

# family classification based on our D/E enrichment findings
DE_ENRICHED = {'RRM_1', 'LSM', 'DEAD', 'Helicase_C'}  # showed D/E enrichment
NON_DE_ENRICHED = {'KH_1', 'S1', 'CSD', 'dsRBD', 'zf-CCHC', 'zf-CCCH'}  # no D/E enrichment or lower

def load_domains():
    """load domain data"""
    with open(RESULTS_DIR / "rbp_domains.json") as f:
        return json.load(f)

def calculate_de_content(seq):
    """calculate D/E percentage"""
    if not seq:
        return 0
    de_count = sum(1 for aa in seq if aa in 'DE')
    return de_count / len(seq) * 100

def find_rgg_motifs(seq):
    """
    find RGG/RG motifs - common in RNA-binding proteins
    patterns: RGG, RG repeats, GRG, etc.
    """
    # RGG repeats
    rgg_pattern = re.compile(r'R[GS]G')
    rgg_matches = rgg_pattern.findall(seq)

    # consecutive RG pairs
    rg_repeats = re.compile(r'(RG){2,}')
    rg_matches = rg_repeats.findall(seq)

    # FG repeats (sometimes in RBPs)
    fg_pattern = re.compile(r'(FG){2,}')
    fg_matches = fg_pattern.findall(seq)

    return {
        'rgg_count': len(rgg_matches),
        'rg_repeat_count': len(rg_matches),
        'fg_repeat_count': len(fg_matches),
        'total_rg_motifs': len(rgg_matches) + len(rg_matches)
    }

def find_sr_motifs(seq):
    """
    find SR-rich regions - common in splicing factors
    patterns: RS, SR repeats
    """
    sr_pattern = re.compile(r'[SR]{3,}')  # 3+ consecutive S or R
    sr_matches = sr_pattern.findall(seq)

    # specific RS dipeptide repeats
    rs_repeats = re.compile(r'(RS){2,}|(SR){2,}')
    rs_matches = rs_repeats.findall(seq)

    return {
        'sr_rich_regions': len(sr_matches),
        'rs_dipeptide_repeats': len(rs_matches),
        'sr_content': (seq.count('S') + seq.count('R')) / len(seq) * 100 if seq else 0
    }

def find_pq_rich(seq):
    """
    find proline-glutamine rich regions - associated with aggregation
    """
    pq_pattern = re.compile(r'[PQ]{3,}')
    pq_matches = pq_pattern.findall(seq)

    return {
        'pq_rich_regions': len(pq_matches),
        'pq_content': (seq.count('P') + seq.count('Q')) / len(seq) * 100 if seq else 0
    }

def calculate_idr_composition(seq):
    """
    calculate IDR-associated amino acid composition
    disorder-promoting: P, E, S, Q, K, A, G
    order-promoting: W, Y, F, I, L, V, C, N
    """
    if not seq:
        return {}

    disorder_promoting = set('PESQKAG')
    order_promoting = set('WYFILVCN')

    disorder_count = sum(1 for aa in seq if aa in disorder_promoting)
    order_count = sum(1 for aa in seq if aa in order_promoting)

    return {
        'disorder_promoting_pct': disorder_count / len(seq) * 100,
        'order_promoting_pct': order_count / len(seq) * 100,
        'disorder_order_ratio': disorder_count / order_count if order_count > 0 else float('inf'),
        # individual important residues
        'glycine_pct': seq.count('G') / len(seq) * 100,
        'proline_pct': seq.count('P') / len(seq) * 100,
        'serine_pct': seq.count('S') / len(seq) * 100,
        'arginine_pct': seq.count('R') / len(seq) * 100,
        'lysine_pct': seq.count('K') / len(seq) * 100,
    }

def analyze_family(domains, family_name):
    """analyze all domains in a family"""
    family_domains = [d for d in domains if d['pfam_name'] == family_name]

    if not family_domains:
        return None

    results = {
        'family': family_name,
        'n_domains': len(family_domains),
        'de_content': [],
        'rgg_motifs': [],
        'sr_motifs': [],
        'pq_motifs': [],
        'idr_composition': []
    }

    for d in family_domains:
        seq = d['sequence']
        results['de_content'].append(calculate_de_content(seq))
        results['rgg_motifs'].append(find_rgg_motifs(seq))
        results['sr_motifs'].append(find_sr_motifs(seq))
        results['pq_motifs'].append(find_pq_rich(seq))
        results['idr_composition'].append(calculate_idr_composition(seq))

    return results

def summarize_family(family_results):
    """summarize family results"""
    if not family_results:
        return None

    n = family_results['n_domains']

    summary = {
        'family': family_results['family'],
        'n_domains': n,
        'mean_de_content': np.mean(family_results['de_content']),
        'std_de_content': np.std(family_results['de_content']),

        # RGG motifs
        'mean_rgg_count': np.mean([m['rgg_count'] for m in family_results['rgg_motifs']]),
        'pct_with_rgg': sum(1 for m in family_results['rgg_motifs'] if m['rgg_count'] > 0) / n * 100,

        # SR motifs
        'mean_sr_content': np.mean([m['sr_content'] for m in family_results['sr_motifs']]),
        'pct_with_sr_regions': sum(1 for m in family_results['sr_motifs'] if m['sr_rich_regions'] > 0) / n * 100,

        # PQ motifs
        'mean_pq_content': np.mean([m['pq_content'] for m in family_results['pq_motifs']]),

        # IDR composition
        'mean_disorder_promoting': np.mean([c['disorder_promoting_pct'] for c in family_results['idr_composition']]),
        'mean_glycine': np.mean([c['glycine_pct'] for c in family_results['idr_composition']]),
        'mean_proline': np.mean([c['proline_pct'] for c in family_results['idr_composition']]),
        'mean_arginine': np.mean([c['arginine_pct'] for c in family_results['idr_composition']]),
    }

    return summary

def compare_groups(de_enriched_data, non_de_data, metric_name):
    """statistical comparison between groups"""
    if not de_enriched_data or not non_de_data:
        return None

    # mann-whitney u test (non-parametric)
    stat, pval = stats.mannwhitneyu(de_enriched_data, non_de_data, alternative='two-sided')

    return {
        'metric': metric_name,
        'de_enriched_mean': np.mean(de_enriched_data),
        'de_enriched_std': np.std(de_enriched_data),
        'non_de_mean': np.mean(non_de_data),
        'non_de_std': np.std(non_de_data),
        'difference': np.mean(de_enriched_data) - np.mean(non_de_data),
        'mann_whitney_p': pval
    }

def main():
    print("=" * 60)
    print("Alternative Regulatory Mechanisms Analysis")
    print("=" * 60)

    # load data
    domains = load_domains()
    print(f"\nLoaded {len(domains)} domains")

    # get available families
    available_families = set(d['pfam_name'] for d in domains)
    print(f"Available families: {len(available_families)}")

    # filter to families we have
    de_families = DE_ENRICHED & available_families
    non_de_families = NON_DE_ENRICHED & available_families

    print(f"\nD/E-enriched families present: {de_families}")
    print(f"Non-D/E families present: {non_de_families}")

    # analyze each family
    all_summaries = []

    print("\n" + "-" * 60)
    print("Per-Family Analysis")
    print("-" * 60)

    for family in sorted(de_families | non_de_families):
        results = analyze_family(domains, family)
        if results:
            summary = summarize_family(results)
            all_summaries.append(summary)

            group = "D/E-enriched" if family in de_families else "Non-D/E"
            print(f"\n{family} ({group}, n={summary['n_domains']}):")
            print(f"  D/E content: {summary['mean_de_content']:.1f}% ± {summary['std_de_content']:.1f}")
            print(f"  RGG motifs: {summary['mean_rgg_count']:.2f} avg, {summary['pct_with_rgg']:.1f}% have any")
            print(f"  SR content: {summary['mean_sr_content']:.1f}%")
            print(f"  Glycine: {summary['mean_glycine']:.1f}%")
            print(f"  Arginine: {summary['mean_arginine']:.1f}%")

    # group comparisons
    print("\n" + "=" * 60)
    print("Group Comparisons: D/E-enriched vs Non-D/E families")
    print("=" * 60)

    # collect all domain-level data for each group
    de_domains = [d for d in domains if d['pfam_name'] in de_families]
    non_de_domains = [d for d in domains if d['pfam_name'] in non_de_families]

    print(f"\nD/E-enriched group: {len(de_domains)} domains from {len(de_families)} families")
    print(f"Non-D/E group: {len(non_de_domains)} domains from {len(non_de_families)} families")

    # compare metrics
    comparisons = []

    # D/E content
    de_de_content = [calculate_de_content(d['sequence']) for d in de_domains]
    non_de_de_content = [calculate_de_content(d['sequence']) for d in non_de_domains]
    comparisons.append(compare_groups(de_de_content, non_de_de_content, 'D/E content (%)'))

    # RGG motifs
    de_rgg = [find_rgg_motifs(d['sequence'])['rgg_count'] for d in de_domains]
    non_de_rgg = [find_rgg_motifs(d['sequence'])['rgg_count'] for d in non_de_domains]
    comparisons.append(compare_groups(de_rgg, non_de_rgg, 'RGG motif count'))

    # glycine content
    de_gly = [d['sequence'].count('G') / len(d['sequence']) * 100 for d in de_domains]
    non_de_gly = [d['sequence'].count('G') / len(d['sequence']) * 100 for d in non_de_domains]
    comparisons.append(compare_groups(de_gly, non_de_gly, 'Glycine (%)'))

    # arginine content
    de_arg = [d['sequence'].count('R') / len(d['sequence']) * 100 for d in de_domains]
    non_de_arg = [d['sequence'].count('R') / len(d['sequence']) * 100 for d in non_de_domains]
    comparisons.append(compare_groups(de_arg, non_de_arg, 'Arginine (%)'))

    # proline content
    de_pro = [d['sequence'].count('P') / len(d['sequence']) * 100 for d in de_domains]
    non_de_pro = [d['sequence'].count('P') / len(d['sequence']) * 100 for d in non_de_domains]
    comparisons.append(compare_groups(de_pro, non_de_pro, 'Proline (%)'))

    # SR content
    de_sr = [find_sr_motifs(d['sequence'])['sr_content'] for d in de_domains]
    non_de_sr = [find_sr_motifs(d['sequence'])['sr_content'] for d in non_de_domains]
    comparisons.append(compare_groups(de_sr, non_de_sr, 'S+R content (%)'))

    # print comparisons
    print("\n" + "-" * 60)
    print(f"{'Metric':<25} {'D/E-enriched':>12} {'Non-D/E':>12} {'Diff':>8} {'p-value':>12}")
    print("-" * 60)

    for comp in comparisons:
        if comp:
            sig = "***" if comp['mann_whitney_p'] < 0.001 else "**" if comp['mann_whitney_p'] < 0.01 else "*" if comp['mann_whitney_p'] < 0.05 else ""
            print(f"{comp['metric']:<25} {comp['de_enriched_mean']:>10.2f} {comp['non_de_mean']:>12.2f} {comp['difference']:>+8.2f} {comp['mann_whitney_p']:>10.2e} {sig}")

    # key findings
    print("\n" + "=" * 60)
    print("KEY FINDINGS")
    print("=" * 60)

    findings = []
    for comp in comparisons:
        if comp and comp['mann_whitney_p'] < 0.05:
            direction = "higher" if comp['difference'] > 0 else "lower"
            findings.append(f"- {comp['metric']}: D/E-enriched families have {direction} values (Δ={comp['difference']:+.2f}, p={comp['mann_whitney_p']:.2e})")

    if findings:
        print("\nSignificant differences:")
        for f in findings:
            print(f)
    else:
        print("\nNo significant differences found at p<0.05")

    # save results
    output = {
        'family_summaries': all_summaries,
        'group_comparisons': comparisons,
        'de_enriched_families': list(de_families),
        'non_de_families': list(non_de_families)
    }

    with open(RESULTS_DIR / "alternative_mechanisms_analysis.json", 'w') as f:
        json.dump(output, f, indent=2, default=str)

    print(f"\nResults saved to {RESULTS_DIR / 'alternative_mechanisms_analysis.json'}")

    return comparisons

if __name__ == "__main__":
    main()
