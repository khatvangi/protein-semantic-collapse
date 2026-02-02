#!/usr/bin/env python3
"""
Three-part analysis of RBD vs non-RBD domain composition:

1. IDR-like region detection (compositional heuristics)
2. D-rich segment mapping (autoinhibitory region hypothesis)
3. Position-specific AA frequency analysis

Tests the hypothesis from Wang et al. 2025 that D/E-rich IDRs
in RBDs serve autoinhibitory functions.
"""

import json
import numpy as np
from collections import Counter, defaultdict
from pathlib import Path
import csv

STANDARD_AAS = "ACDEFGHIKLMNPQRSTVWY"
STANDARD_AAS_SET = set(STANDARD_AAS)

# charged residues
ACIDIC = set("DE")
BASIC = set("KR")
CHARGED = ACIDIC | BASIC

# disorder-promoting residues (from literature)
DISORDER_PROMOTING = set("AEGRQSPK")
ORDER_PROMOTING = set("WFYILMVNC")


def compute_local_entropy(window):
    """Shannon entropy of a short window (for low-complexity detection)."""
    if len(window) < 3:
        return 0.0
    counts = Counter(window)
    total = len(window)
    entropy = 0.0
    for c in counts.values():
        p = c / total
        if p > 0:
            entropy -= p * np.log2(p)
    return entropy


def find_de_rich_segments(sequence, window_size=10, threshold=0.4):
    """
    Find D/E-rich segments using sliding window.

    Args:
        sequence: amino acid sequence
        window_size: sliding window size
        threshold: minimum fraction of D/E to count as enriched

    Returns:
        list of (start, end, de_fraction) tuples
    """
    segments = []
    if len(sequence) < window_size:
        return segments

    in_segment = False
    seg_start = 0

    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i + window_size]
        de_count = sum(1 for aa in window if aa in ACIDIC)
        de_frac = de_count / window_size

        if de_frac >= threshold:
            if not in_segment:
                seg_start = i
                in_segment = True
        else:
            if in_segment:
                # end of segment - calculate overall D/E content
                seg_seq = sequence[seg_start:i + window_size - 1]
                total_de = sum(1 for aa in seg_seq if aa in ACIDIC)
                segments.append((seg_start, i + window_size - 1, total_de / len(seg_seq)))
                in_segment = False

    # handle segment at end
    if in_segment:
        seg_seq = sequence[seg_start:]
        total_de = sum(1 for aa in seg_seq if aa in ACIDIC)
        segments.append((seg_start, len(sequence), total_de / len(seg_seq)))

    return segments


def compute_disorder_score(sequence, window_size=15):
    """
    Simple compositional disorder score.

    Returns score between -1 (ordered) and +1 (disordered).
    Based on disorder-promoting vs order-promoting residue balance.
    """
    if len(sequence) < window_size:
        return 0.0

    scores = []
    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i + window_size]
        n_disorder = sum(1 for aa in window if aa in DISORDER_PROMOTING)
        n_order = sum(1 for aa in window if aa in ORDER_PROMOTING)
        # normalize
        score = (n_disorder - n_order) / window_size
        scores.append(score)

    return np.mean(scores) if scores else 0.0


def analyze_position_specific(sequence, n_bins=5):
    """
    Compute AA frequencies by normalized position within domain.

    Bins: N-term (0-0.2), early-mid (0.2-0.4), middle (0.4-0.6),
          late-mid (0.6-0.8), C-term (0.8-1.0)

    Returns dict {bin_idx: {AA: count}}
    """
    bin_counts = {i: Counter() for i in range(n_bins)}

    n = len(sequence)
    if n < n_bins:
        return bin_counts

    for pos, aa in enumerate(sequence):
        if aa not in STANDARD_AAS_SET:
            continue
        # normalized position [0, 1)
        norm_pos = pos / n
        bin_idx = min(int(norm_pos * n_bins), n_bins - 1)
        bin_counts[bin_idx][aa] += 1

    return bin_counts


def main():
    base_dir = Path("/storage/kiran-stuff/protein-semantic-collapse")
    data_dir = base_dir / "data"
    out_dir = base_dir / "results"
    out_dir.mkdir(exist_ok=True)

    domain_file = data_dir / "domain_sequences.jsonl"

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

    print(f"  RBD domains: {len(rbd_domains):,}")
    print(f"  non-RBD domains: {len(nonrbd_domains):,}")

    # =========================================================
    # ANALYSIS 1: D/E-RICH SEGMENT DETECTION
    # =========================================================
    print("\n" + "=" * 70)
    print("ANALYSIS 1: D/E-RICH SEGMENT MAPPING")
    print("=" * 70)

    def analyze_de_segments(domains, label):
        """Analyze D/E-rich segments in domain set."""
        results = {
            "n_domains": len(domains),
            "n_with_de_segment": 0,
            "total_segments": 0,
            "segment_lengths": [],
            "segment_densities": [],  # D/E fraction within segment
            "segment_positions": [],  # normalized start position
            "domains_with_segments": []
        }

        for d in domains:
            seq = d["sequence"]
            segments = find_de_rich_segments(seq, window_size=8, threshold=0.375)  # 3/8 D/E

            if segments:
                results["n_with_de_segment"] += 1
                results["total_segments"] += len(segments)

                for start, end, density in segments:
                    results["segment_lengths"].append(end - start)
                    results["segment_densities"].append(density)
                    results["segment_positions"].append(start / len(seq))  # normalized

                results["domains_with_segments"].append({
                    "accession": d["accession"],
                    "pfam": d["pfam"],
                    "n_segments": len(segments),
                    "total_de_length": sum(e - s for s, e, _ in segments)
                })

        return results

    rbd_de = analyze_de_segments(rbd_domains, "RBD")
    nonrbd_de = analyze_de_segments(nonrbd_domains, "non-RBD")

    print(f"\n{'Metric':<40} {'RBD':>15} {'non-RBD':>15}")
    print("-" * 70)
    print(f"{'Total domains':<40} {rbd_de['n_domains']:>15,} {nonrbd_de['n_domains']:>15,}")
    print(f"{'Domains with D/E-rich segment':<40} {rbd_de['n_with_de_segment']:>15,} {nonrbd_de['n_with_de_segment']:>15,}")

    rbd_pct = 100 * rbd_de['n_with_de_segment'] / rbd_de['n_domains']
    nonrbd_pct = 100 * nonrbd_de['n_with_de_segment'] / nonrbd_de['n_domains']
    print(f"{'% with D/E segment':<40} {rbd_pct:>14.1f}% {nonrbd_pct:>14.1f}%")

    if rbd_de['segment_lengths']:
        print(f"{'Mean segment length (residues)':<40} {np.mean(rbd_de['segment_lengths']):>15.1f} {np.mean(nonrbd_de['segment_lengths']):>15.1f}")
        print(f"{'Mean D/E density in segments':<40} {np.mean(rbd_de['segment_densities']):>15.3f} {np.mean(nonrbd_de['segment_densities']):>15.3f}")

    # position distribution of segments
    if rbd_de['segment_positions']:
        print(f"\nSegment position distribution (normalized 0=N-term, 1=C-term):")
        rbd_pos = np.array(rbd_de['segment_positions'])
        nonrbd_pos = np.array(nonrbd_de['segment_positions']) if nonrbd_de['segment_positions'] else np.array([0.5])
        print(f"  RBD:     mean={np.mean(rbd_pos):.3f}, median={np.median(rbd_pos):.3f}")
        print(f"  non-RBD: mean={np.mean(nonrbd_pos):.3f}, median={np.median(nonrbd_pos):.3f}")

        # terminal enrichment
        rbd_nterm = np.sum(rbd_pos < 0.2) / len(rbd_pos) * 100
        rbd_cterm = np.sum(rbd_pos > 0.8) / len(rbd_pos) * 100
        print(f"  RBD N-terminal (<0.2): {rbd_nterm:.1f}%")
        print(f"  RBD C-terminal (>0.8): {rbd_cterm:.1f}%")

    # =========================================================
    # ANALYSIS 2: DISORDER SCORE COMPARISON
    # =========================================================
    print("\n" + "=" * 70)
    print("ANALYSIS 2: COMPOSITIONAL DISORDER SCORES")
    print("=" * 70)

    rbd_disorder = [compute_disorder_score(d["sequence"]) for d in rbd_domains]
    nonrbd_disorder = [compute_disorder_score(d["sequence"]) for d in nonrbd_domains]

    print(f"\n{'Metric':<40} {'RBD':>15} {'non-RBD':>15}")
    print("-" * 70)
    print(f"{'Mean disorder score':<40} {np.mean(rbd_disorder):>15.4f} {np.mean(nonrbd_disorder):>15.4f}")
    print(f"{'Median disorder score':<40} {np.median(rbd_disorder):>15.4f} {np.median(nonrbd_disorder):>15.4f}")
    print(f"{'Std disorder score':<40} {np.std(rbd_disorder):>15.4f} {np.std(nonrbd_disorder):>15.4f}")

    # fraction with positive (disordered) score
    rbd_disordered = sum(1 for s in rbd_disorder if s > 0) / len(rbd_disorder) * 100
    nonrbd_disordered = sum(1 for s in nonrbd_disorder if s > 0) / len(nonrbd_disorder) * 100
    print(f"{'% with disorder bias (score > 0)':<40} {rbd_disordered:>14.1f}% {nonrbd_disordered:>14.1f}%")

    # =========================================================
    # ANALYSIS 3: POSITION-SPECIFIC AA FREQUENCIES
    # =========================================================
    print("\n" + "=" * 70)
    print("ANALYSIS 3: POSITION-SPECIFIC AA FREQUENCIES")
    print("=" * 70)

    N_BINS = 5
    BIN_NAMES = ["N-term", "Early-mid", "Middle", "Late-mid", "C-term"]

    # accumulate counts
    rbd_pos_counts = {i: Counter() for i in range(N_BINS)}
    nonrbd_pos_counts = {i: Counter() for i in range(N_BINS)}

    for d in rbd_domains:
        bin_counts = analyze_position_specific(d["sequence"], N_BINS)
        for i, counts in bin_counts.items():
            rbd_pos_counts[i].update(counts)

    for d in nonrbd_domains:
        bin_counts = analyze_position_specific(d["sequence"], N_BINS)
        for i, counts in bin_counts.items():
            nonrbd_pos_counts[i].update(counts)

    # convert to frequencies
    def counts_to_freqs(counts_dict):
        freqs = {}
        for bin_idx, counts in counts_dict.items():
            total = sum(counts.values())
            if total > 0:
                freqs[bin_idx] = {aa: counts.get(aa, 0) / total for aa in STANDARD_AAS}
            else:
                freqs[bin_idx] = {aa: 0.0 for aa in STANDARD_AAS}
        return freqs

    rbd_pos_freqs = counts_to_freqs(rbd_pos_counts)
    nonrbd_pos_freqs = counts_to_freqs(nonrbd_pos_counts)

    # focus on key residues: D, E, K, R (charged), F, Y (aromatic)
    KEY_AAS = ["D", "E", "K", "R", "F", "Y"]

    print(f"\nKey AA frequencies by position (%):")
    print(f"{'Position':<12}", end="")
    for aa in KEY_AAS:
        print(f"  {aa}-RBD  {aa}-non", end="")
    print()
    print("-" * 90)

    for i in range(N_BINS):
        print(f"{BIN_NAMES[i]:<12}", end="")
        for aa in KEY_AAS:
            rbd_f = rbd_pos_freqs[i][aa] * 100
            nonrbd_f = nonrbd_pos_freqs[i][aa] * 100
            print(f"  {rbd_f:5.2f}  {nonrbd_f:5.2f}", end="")
        print()

    # D+E combined (acidic)
    print(f"\n{'Acidic (D+E) by position:'}")
    print(f"{'Position':<12} {'RBD %':>10} {'non-RBD %':>12} {'Diff':>10}")
    print("-" * 50)
    for i in range(N_BINS):
        rbd_de = (rbd_pos_freqs[i]["D"] + rbd_pos_freqs[i]["E"]) * 100
        nonrbd_de = (nonrbd_pos_freqs[i]["D"] + nonrbd_pos_freqs[i]["E"]) * 100
        diff = rbd_de - nonrbd_de
        marker = "**" if abs(diff) > 0.5 else ""
        print(f"{BIN_NAMES[i]:<12} {rbd_de:>10.2f} {nonrbd_de:>12.2f} {diff:>+10.2f} {marker}")

    # K+R combined (basic)
    print(f"\n{'Basic (K+R) by position:'}")
    print(f"{'Position':<12} {'RBD %':>10} {'non-RBD %':>12} {'Diff':>10}")
    print("-" * 50)
    for i in range(N_BINS):
        rbd_kr = (rbd_pos_freqs[i]["K"] + rbd_pos_freqs[i]["R"]) * 100
        nonrbd_kr = (nonrbd_pos_freqs[i]["K"] + nonrbd_pos_freqs[i]["R"]) * 100
        diff = rbd_kr - nonrbd_kr
        marker = "**" if abs(diff) > 0.5 else ""
        print(f"{BIN_NAMES[i]:<12} {rbd_kr:>10.2f} {nonrbd_kr:>12.2f} {diff:>+10.2f} {marker}")

    # net charge by position
    print(f"\n{'Net charge (K+R-D-E) by position:'}")
    print(f"{'Position':<12} {'RBD':>10} {'non-RBD':>12} {'Diff':>10}")
    print("-" * 50)
    for i in range(N_BINS):
        rbd_net = (rbd_pos_freqs[i]["K"] + rbd_pos_freqs[i]["R"] -
                   rbd_pos_freqs[i]["D"] - rbd_pos_freqs[i]["E"]) * 100
        nonrbd_net = (nonrbd_pos_freqs[i]["K"] + nonrbd_pos_freqs[i]["R"] -
                      nonrbd_pos_freqs[i]["D"] - nonrbd_pos_freqs[i]["E"]) * 100
        diff = rbd_net - nonrbd_net
        print(f"{BIN_NAMES[i]:<12} {rbd_net:>+10.2f} {nonrbd_net:>+12.2f} {diff:>+10.2f}")

    # =========================================================
    # SAVE DETAILED RESULTS
    # =========================================================

    # D/E segment details
    with open(out_dir / "de_rich_segments.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["domain_type", "accession", "pfam", "n_segments", "total_de_length"])
        for d in rbd_de["domains_with_segments"]:
            writer.writerow(["RBD", d["accession"], d["pfam"], d["n_segments"], d["total_de_length"]])
        for d in nonrbd_de["domains_with_segments"]:
            writer.writerow(["non-RBD", d["accession"], d["pfam"], d["n_segments"], d["total_de_length"]])

    # position-specific frequencies
    with open(out_dir / "position_specific_freqs.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["domain_type", "position", "AA", "frequency"])
        for i in range(N_BINS):
            for aa in STANDARD_AAS:
                writer.writerow(["RBD", BIN_NAMES[i], aa, f"{rbd_pos_freqs[i][aa]:.6f}"])
                writer.writerow(["non-RBD", BIN_NAMES[i], aa, f"{nonrbd_pos_freqs[i][aa]:.6f}"])

    # disorder scores
    with open(out_dir / "disorder_scores.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["domain_type", "accession", "pfam", "disorder_score"])
        for d, score in zip(rbd_domains, rbd_disorder):
            writer.writerow(["RBD", d["accession"], d["pfam"], f"{score:.4f}"])
        for d, score in zip(nonrbd_domains, nonrbd_disorder):
            writer.writerow(["non-RBD", d["accession"], d["pfam"], f"{score:.4f}"])

    print(f"\n\nResults saved to {out_dir}/")
    print("  - de_rich_segments.csv")
    print("  - position_specific_freqs.csv")
    print("  - disorder_scores.csv")

    # =========================================================
    # SUMMARY: TEST AUTOINHIBITION HYPOTHESIS
    # =========================================================
    print("\n" + "=" * 70)
    print("SUMMARY: AUTOINHIBITION HYPOTHESIS TEST")
    print("=" * 70)

    print("\nHypothesis (Wang et al. 2025): RBDs contain D/E-rich IDRs that")
    print("function as autoinhibitory elements, competing with RNA for binding.")

    print("\nEvidence from our data:")

    # test 1: more D/E segments in RBDs?
    if rbd_pct > nonrbd_pct:
        print(f"  [+] RBDs have MORE D/E-rich segments ({rbd_pct:.1f}% vs {nonrbd_pct:.1f}%)")
    else:
        print(f"  [-] RBDs have FEWER D/E-rich segments ({rbd_pct:.1f}% vs {nonrbd_pct:.1f}%)")

    # test 2: higher disorder scores in RBDs?
    if np.mean(rbd_disorder) > np.mean(nonrbd_disorder):
        print(f"  [+] RBDs have HIGHER disorder scores ({np.mean(rbd_disorder):.4f} vs {np.mean(nonrbd_disorder):.4f})")
    else:
        print(f"  [-] RBDs have LOWER disorder scores ({np.mean(rbd_disorder):.4f} vs {np.mean(nonrbd_disorder):.4f})")

    # test 3: D/E enrichment at termini?
    if rbd_de['segment_positions']:
        terminal_enrichment = (rbd_nterm + rbd_cterm)
        print(f"  [?] D/E segments at termini: {terminal_enrichment:.1f}% (N={rbd_nterm:.1f}%, C={rbd_cterm:.1f}%)")

    # test 4: more acidic at specific positions?
    max_diff_pos = max(range(N_BINS), key=lambda i:
        (rbd_pos_freqs[i]["D"] + rbd_pos_freqs[i]["E"]) -
        (nonrbd_pos_freqs[i]["D"] + nonrbd_pos_freqs[i]["E"]))
    max_diff = ((rbd_pos_freqs[max_diff_pos]["D"] + rbd_pos_freqs[max_diff_pos]["E"]) -
                (nonrbd_pos_freqs[max_diff_pos]["D"] + nonrbd_pos_freqs[max_diff_pos]["E"])) * 100
    print(f"  [?] Largest D+E enrichment in RBDs: {BIN_NAMES[max_diff_pos]} ({max_diff:+.2f}%)")


if __name__ == "__main__":
    main()
