#!/usr/bin/env python3
"""
Analyze amino acid vocabulary usage in RBPs vs non-RBPs.

Tests the hypothesis that RBPs use a reduced AA vocabulary
(semantic collapse) compared to general proteins.
"""

import numpy as np
from collections import Counter
from pathlib import Path
import csv
from scipy import stats

# standard 20 amino acids
STANDARD_AAS = set("ACDEFGHIKLMNPQRSTVWY")

def parse_fasta(fasta_path):
    """yield (id, sequence) tuples from fasta file"""
    current_id = None
    current_seq = []

    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id:
                    yield current_id, "".join(current_seq)
                # extract uniprot id (first part after >)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.upper())

        if current_id:
            yield current_id, "".join(current_seq)

def compute_aa_frequencies(sequence):
    """compute frequencies of standard AAs in sequence"""
    # filter to standard AAs only
    filtered = [aa for aa in sequence if aa in STANDARD_AAS]
    if len(filtered) == 0:
        return None

    counts = Counter(filtered)
    total = sum(counts.values())

    # return frequencies for all 20 AAs (0 if not present)
    freqs = {aa: counts.get(aa, 0) / total for aa in STANDARD_AAS}
    return freqs

def compute_entropy(freq_dict):
    """compute Shannon entropy from frequency dict"""
    # H = -sum(p * log2(p)) for p > 0
    entropy = 0
    for p in freq_dict.values():
        if p > 0:
            entropy -= p * np.log2(p)
    return entropy

def effective_alphabet_size(entropy):
    """convert entropy to effective alphabet size: 2^H"""
    return 2 ** entropy

def analyze_dataset(fasta_path, label):
    """analyze all proteins in a fasta file"""
    results = []
    all_aa_counts = Counter()

    for prot_id, seq in parse_fasta(fasta_path):
        freqs = compute_aa_frequencies(seq)
        if freqs is None:
            continue

        # per-protein entropy
        ent = compute_entropy(freqs)
        eff_size = effective_alphabet_size(ent)

        # count AAs for aggregate frequencies
        for aa in seq:
            if aa in STANDARD_AAS:
                all_aa_counts[aa] += 1

        results.append({
            "id": prot_id,
            "label": label,
            "length": len([aa for aa in seq if aa in STANDARD_AAS]),
            "entropy": ent,
            "effective_alphabet": eff_size,
            **{f"freq_{aa}": freqs[aa] for aa in sorted(STANDARD_AAS)}
        })

    # aggregate frequencies
    total = sum(all_aa_counts.values())
    agg_freqs = {aa: all_aa_counts.get(aa, 0) / total for aa in STANDARD_AAS}
    agg_entropy = compute_entropy(agg_freqs)

    return results, agg_freqs, agg_entropy

def main():
    base_dir = Path("/storage/kiran-stuff/protein-semantic-collapse")
    seq_dir = base_dir / "sequences"
    out_dir = base_dir / "results"
    out_dir.mkdir(exist_ok=True)

    print("Analyzing RBPs...")
    rbp_results, rbp_agg_freqs, rbp_agg_entropy = analyze_dataset(
        seq_dir / "rbp_sequences.fasta", "RBP"
    )
    print(f"  {len(rbp_results)} proteins analyzed")

    print("Analyzing non-RBPs...")
    nonrbp_results, nonrbp_agg_freqs, nonrbp_agg_entropy = analyze_dataset(
        seq_dir / "nonrbp_sequences.fasta", "non-RBP"
    )
    print(f"  {len(nonrbp_results)} proteins analyzed")

    # save per-protein results
    all_results = rbp_results + nonrbp_results
    fieldnames = list(all_results[0].keys())

    with open(out_dir / "per_protein_stats.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_results)
    print(f"\nPer-protein stats saved to {out_dir / 'per_protein_stats.csv'}")

    # extract entropy distributions
    rbp_entropies = [r["entropy"] for r in rbp_results]
    nonrbp_entropies = [r["entropy"] for r in nonrbp_results]

    # statistical comparison
    stat, pvalue = stats.mannwhitneyu(rbp_entropies, nonrbp_entropies, alternative="two-sided")

    # effect size (rank-biserial correlation)
    n1, n2 = len(rbp_entropies), len(nonrbp_entropies)
    effect_size = 1 - (2 * stat) / (n1 * n2)

    # summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)

    print(f"\n{'Metric':<30} {'RBP':>12} {'non-RBP':>12}")
    print("-" * 54)
    print(f"{'N proteins':<30} {len(rbp_results):>12,} {len(nonrbp_results):>12,}")
    print(f"{'Mean entropy (bits)':<30} {np.mean(rbp_entropies):>12.4f} {np.mean(nonrbp_entropies):>12.4f}")
    print(f"{'Median entropy (bits)':<30} {np.median(rbp_entropies):>12.4f} {np.median(nonrbp_entropies):>12.4f}")
    print(f"{'Std entropy':<30} {np.std(rbp_entropies):>12.4f} {np.std(nonrbp_entropies):>12.4f}")
    print(f"{'Aggregate entropy (bits)':<30} {rbp_agg_entropy:>12.4f} {nonrbp_agg_entropy:>12.4f}")
    print(f"{'Effective alphabet (agg)':<30} {effective_alphabet_size(rbp_agg_entropy):>12.2f} {effective_alphabet_size(nonrbp_agg_entropy):>12.2f}")

    print(f"\n{'Statistical comparison':}")
    print(f"  Mann-Whitney U: {stat:,.0f}")
    print(f"  p-value: {pvalue:.2e}")
    print(f"  Effect size (rank-biserial r): {effect_size:.4f}")

    # interpretation
    print("\n" + "=" * 60)
    print("INTERPRETATION")
    print("=" * 60)
    entropy_diff = np.mean(rbp_entropies) - np.mean(nonrbp_entropies)
    if entropy_diff < 0:
        print(f"\nRBPs have LOWER entropy by {abs(entropy_diff):.4f} bits")
        print("→ Supports semantic collapse hypothesis")
    else:
        print(f"\nRBPs have HIGHER entropy by {entropy_diff:.4f} bits")
        print("→ Does NOT support semantic collapse hypothesis")

    if abs(effect_size) < 0.1:
        print(f"Effect size is negligible (r = {effect_size:.4f})")
    elif abs(effect_size) < 0.3:
        print(f"Effect size is small (r = {effect_size:.4f})")
    elif abs(effect_size) < 0.5:
        print(f"Effect size is medium (r = {effect_size:.4f})")
    else:
        print(f"Effect size is large (r = {effect_size:.4f})")

    # save aggregate frequencies for plotting
    with open(out_dir / "aggregate_frequencies.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["AA", "RBP_freq", "nonRBP_freq", "diff"])
        for aa in sorted(STANDARD_AAS):
            diff = rbp_agg_freqs[aa] - nonrbp_agg_freqs[aa]
            writer.writerow([aa, f"{rbp_agg_freqs[aa]:.6f}", f"{nonrbp_agg_freqs[aa]:.6f}", f"{diff:.6f}"])
    print(f"\nAggregate frequencies saved to {out_dir / 'aggregate_frequencies.csv'}")

    # top AA differences
    print("\n" + "=" * 60)
    print("TOP AA FREQUENCY DIFFERENCES (RBP - non-RBP)")
    print("=" * 60)
    diffs = [(aa, rbp_agg_freqs[aa] - nonrbp_agg_freqs[aa]) for aa in STANDARD_AAS]
    diffs.sort(key=lambda x: x[1])

    print("\nMost depleted in RBPs:")
    for aa, d in diffs[:5]:
        print(f"  {aa}: {d:+.4f} ({rbp_agg_freqs[aa]:.4f} vs {nonrbp_agg_freqs[aa]:.4f})")

    print("\nMost enriched in RBPs:")
    for aa, d in diffs[-5:][::-1]:
        print(f"  {aa}: {d:+.4f} ({rbp_agg_freqs[aa]:.4f} vs {nonrbp_agg_freqs[aa]:.4f})")

if __name__ == "__main__":
    main()
