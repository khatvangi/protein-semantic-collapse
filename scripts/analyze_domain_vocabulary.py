#!/usr/bin/env python3
"""
Analyze AA vocabulary at the domain level.

Compares RBD (RNA-binding domains) vs non-RBD functional domains.
Tests whether RBDs show reduced AA vocabulary (semantic collapse).

Key improvements over protein-level analysis:
1. Uses actual domain sequences (not full proteins)
2. Fair comparison: domain vs domain (not domain vs protein)
3. Controls for domain length effects
4. Stratifies by domain family for deeper insights
"""

import json
import numpy as np
from collections import Counter
from pathlib import Path
import csv
from scipy import stats

STANDARD_AAS = "ACDEFGHIKLMNPQRSTVWY"
STANDARD_AAS_SET = set(STANDARD_AAS)

def compute_aa_frequencies(sequence):
    """
    Compute frequencies of standard 20 AAs.
    Returns dict {AA: freq} or None if sequence too short.
    """
    filtered = [aa for aa in sequence if aa in STANDARD_AAS_SET]
    if len(filtered) < 20:  # need enough residues for meaningful stats
        return None

    counts = Counter(filtered)
    total = sum(counts.values())
    return {aa: counts.get(aa, 0) / total for aa in STANDARD_AAS}

def compute_entropy(freq_dict):
    """Shannon entropy in bits: H = -sum(p * log2(p))"""
    entropy = 0.0
    for p in freq_dict.values():
        if p > 0:
            entropy -= p * np.log2(p)
    return entropy

def effective_alphabet_size(entropy):
    """Convert entropy to effective alphabet size: 2^H"""
    return 2 ** entropy

def analyze_domains(domains):
    """
    Analyze a list of domain dicts.
    Returns (per_domain_results, aggregate_freqs, aggregate_entropy).
    """
    results = []
    all_aa_counts = Counter()
    lengths = []

    for dom in domains:
        seq = dom["sequence"]
        freqs = compute_aa_frequencies(seq)
        if freqs is None:
            continue

        ent = compute_entropy(freqs)
        eff = effective_alphabet_size(ent)

        # accumulate for aggregate stats
        for aa in seq:
            if aa in STANDARD_AAS_SET:
                all_aa_counts[aa] += 1
        lengths.append(dom["domain_length"])

        results.append({
            "accession": dom["accession"],
            "pfam": dom["pfam"],
            "pfam_name": dom["pfam_name"],
            "length": dom["domain_length"],
            "entropy": ent,
            "effective_alphabet": eff,
            "is_rbd": dom["is_rbd"]
        })

    # aggregate frequencies
    total = sum(all_aa_counts.values())
    if total > 0:
        agg_freqs = {aa: all_aa_counts.get(aa, 0) / total for aa in STANDARD_AAS}
        agg_entropy = compute_entropy(agg_freqs)
    else:
        agg_freqs = {aa: 0.0 for aa in STANDARD_AAS}
        agg_entropy = 0.0

    return results, agg_freqs, agg_entropy, lengths

def main():
    base_dir = Path("/storage/kiran-stuff/protein-semantic-collapse")
    data_dir = base_dir / "data"
    out_dir = base_dir / "results"
    out_dir.mkdir(exist_ok=True)

    # load domain sequences
    domain_file = data_dir / "domain_sequences.jsonl"
    if not domain_file.exists():
        print(f"ERROR: {domain_file} not found")
        print("Run fetch_domain_sequences.py first")
        return

    print("Loading domain sequences...")
    rbd_domains = []
    nonrbd_domains = []

    with open(domain_file) as f:
        for line in f:
            d = json.loads(line)
            if d["is_rbd"]:
                rbd_domains.append(d)
            else:
                nonrbd_domains.append(d)

    print(f"  RBD domains: {len(rbd_domains)}")
    print(f"  non-RBD domains: {len(nonrbd_domains)}")

    if not rbd_domains or not nonrbd_domains:
        print("ERROR: Not enough domain data")
        return

    # analyze each class
    print("\nAnalyzing RBD domains...")
    rbd_results, rbd_freqs, rbd_agg_ent, rbd_lengths = analyze_domains(rbd_domains)
    print(f"  Valid domains: {len(rbd_results)}")

    print("Analyzing non-RBD domains...")
    nonrbd_results, nonrbd_freqs, nonrbd_agg_ent, nonrbd_lengths = analyze_domains(nonrbd_domains)
    print(f"  Valid domains: {len(nonrbd_results)}")

    # extract entropies
    rbd_entropies = [r["entropy"] for r in rbd_results]
    nonrbd_entropies = [r["entropy"] for r in nonrbd_results]

    # =========================================================
    # STATISTICAL TESTS
    # =========================================================
    # Mann-Whitney U (non-parametric)
    stat, pvalue = stats.mannwhitneyu(rbd_entropies, nonrbd_entropies, alternative="two-sided")
    n1, n2 = len(rbd_entropies), len(nonrbd_entropies)
    # rank-biserial correlation as effect size
    effect_size = 1 - (2 * stat) / (n1 * n2)

    # Welch's t-test (parametric, for comparison)
    t_stat, t_pvalue = stats.ttest_ind(rbd_entropies, nonrbd_entropies, equal_var=False)

    # Cohen's d effect size
    pooled_std = np.sqrt((np.std(rbd_entropies)**2 + np.std(nonrbd_entropies)**2) / 2)
    cohens_d = (np.mean(rbd_entropies) - np.mean(nonrbd_entropies)) / pooled_std if pooled_std > 0 else 0

    # =========================================================
    # OUTPUT
    # =========================================================
    print("\n" + "=" * 70)
    print("DOMAIN-LEVEL AA VOCABULARY ANALYSIS")
    print("=" * 70)

    print(f"\n{'Metric':<35} {'RBD':>15} {'non-RBD':>15}")
    print("-" * 65)
    print(f"{'N domains analyzed':<35} {len(rbd_results):>15,} {len(nonrbd_results):>15,}")
    print(f"{'Mean length (residues)':<35} {np.mean(rbd_lengths):>15.1f} {np.mean(nonrbd_lengths):>15.1f}")
    print(f"{'Mean entropy (bits)':<35} {np.mean(rbd_entropies):>15.4f} {np.mean(nonrbd_entropies):>15.4f}")
    print(f"{'Median entropy (bits)':<35} {np.median(rbd_entropies):>15.4f} {np.median(nonrbd_entropies):>15.4f}")
    print(f"{'Std entropy':<35} {np.std(rbd_entropies):>15.4f} {np.std(nonrbd_entropies):>15.4f}")
    print(f"{'Aggregate entropy (pooled)':<35} {rbd_agg_ent:>15.4f} {nonrbd_agg_ent:>15.4f}")
    print(f"{'Effective alphabet size (2^H)':<35} {2**rbd_agg_ent:>15.2f} {2**nonrbd_agg_ent:>15.2f}")

    print(f"\n{'STATISTICAL TESTS':}")
    print(f"  Mann-Whitney U: {stat:,.0f}")
    print(f"  p-value (U): {pvalue:.2e}")
    print(f"  Effect size (rank-biserial r): {effect_size:.4f}")
    print(f"  Welch's t: {t_stat:.4f}")
    print(f"  p-value (t): {t_pvalue:.2e}")
    print(f"  Cohen's d: {cohens_d:.4f}")

    # interpretation
    print("\n" + "=" * 70)
    print("INTERPRETATION")
    print("=" * 70)

    diff = np.mean(rbd_entropies) - np.mean(nonrbd_entropies)
    if diff < 0:
        direction = "LOWER"
        hypothesis = "SUPPORTS"
    else:
        direction = "HIGHER"
        hypothesis = "Does NOT support"

    print(f"\nRBDs have {direction} entropy by {abs(diff):.4f} bits (mean)")
    print(f"→ {hypothesis} semantic collapse hypothesis")

    # effect size interpretation
    if abs(cohens_d) < 0.2:
        size_desc = "negligible"
    elif abs(cohens_d) < 0.5:
        size_desc = "small"
    elif abs(cohens_d) < 0.8:
        size_desc = "medium"
    else:
        size_desc = "large"
    print(f"Effect size is {size_desc} (Cohen's d = {cohens_d:.4f})")

    # =========================================================
    # AA FREQUENCY DIFFERENCES
    # =========================================================
    print("\n" + "=" * 70)
    print("TOP AA FREQUENCY DIFFERENCES (RBD - non-RBD)")
    print("=" * 70)

    diffs = [(aa, rbd_freqs[aa] - nonrbd_freqs[aa]) for aa in STANDARD_AAS]
    diffs.sort(key=lambda x: x[1])

    print("\nMost DEPLETED in RBDs:")
    for aa, d in diffs[:5]:
        pct = d * 100
        print(f"  {aa}: {pct:+.2f}% ({rbd_freqs[aa]*100:.2f}% vs {nonrbd_freqs[aa]*100:.2f}%)")

    print("\nMost ENRICHED in RBDs:")
    for aa, d in diffs[-5:][::-1]:
        pct = d * 100
        print(f"  {aa}: {pct:+.2f}% ({rbd_freqs[aa]*100:.2f}% vs {nonrbd_freqs[aa]*100:.2f}%)")

    # =========================================================
    # BY-FAMILY BREAKDOWN
    # =========================================================
    print("\n" + "=" * 70)
    print("ENTROPY BY DOMAIN FAMILY (top 10)")
    print("=" * 70)

    # aggregate by Pfam family
    family_entropies = {}
    for r in rbd_results + nonrbd_results:
        key = (r["pfam"], r["pfam_name"], r["is_rbd"])
        if key not in family_entropies:
            family_entropies[key] = []
        family_entropies[key].append(r["entropy"])

    # summarize
    family_summary = []
    for (pfam, name, is_rbd), entropies in family_entropies.items():
        if len(entropies) >= 10:  # minimum sample size
            family_summary.append({
                "pfam": pfam,
                "name": name,
                "is_rbd": is_rbd,
                "n": len(entropies),
                "mean_entropy": np.mean(entropies),
                "std_entropy": np.std(entropies)
            })

    # sort by mean entropy
    family_summary.sort(key=lambda x: x["mean_entropy"])

    print(f"\n{'Family':<12} {'Name':<20} {'Type':<8} {'N':>6} {'Mean H':>8} {'Std':>6}")
    print("-" * 70)
    for f in family_summary[:10]:
        t = "RBD" if f["is_rbd"] else "non-RBD"
        print(f"{f['pfam']:<12} {f['name'][:20]:<20} {t:<8} {f['n']:>6} {f['mean_entropy']:>8.4f} {f['std_entropy']:>6.4f}")

    # =========================================================
    # SAVE RESULTS
    # =========================================================
    # per-domain stats
    all_results = rbd_results + nonrbd_results
    with open(out_dir / "per_domain_stats.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=all_results[0].keys())
        writer.writeheader()
        writer.writerows(all_results)

    # aggregate frequencies
    with open(out_dir / "domain_aggregate_frequencies.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["AA", "RBD_freq", "nonRBD_freq", "diff", "diff_pct"])
        for aa in STANDARD_AAS:
            d = rbd_freqs[aa] - nonrbd_freqs[aa]
            writer.writerow([aa, f"{rbd_freqs[aa]:.6f}", f"{nonrbd_freqs[aa]:.6f}",
                           f"{d:.6f}", f"{d*100:.4f}"])

    # family summary
    with open(out_dir / "family_entropy_summary.csv", "w", newline="") as f:
        fieldnames = ["pfam", "name", "is_rbd", "n", "mean_entropy", "std_entropy"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(family_summary)

    print(f"\nResults saved to {out_dir}/")
    print("  - per_domain_stats.csv")
    print("  - domain_aggregate_frequencies.csv")
    print("  - family_entropy_summary.csv")

if __name__ == "__main__":
    main()
