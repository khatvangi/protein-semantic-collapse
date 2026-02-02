#!/usr/bin/env python3
"""
Three analyses to find novel insights:
1. Family outliers - RBDs that don't follow D/E pattern
2. Entropy vs autoinhibition correlation
3. Prepare data for structural validation
"""

import json
import numpy as np
from collections import defaultdict
from pathlib import Path
from scipy import stats
from Bio import SeqIO

STANDARD_AAS = set("ACDEFGHIKLMNPQRSTVWY")
ACIDIC = set("DE")
BASIC = set("KR")
DISORDER_PROMOTING = set("AEGRQSPK")
ORDER_PROMOTING = set("WFYILMVNC")


def compute_de_fraction(seq):
    valid = [aa for aa in seq if aa in STANDARD_AAS]
    if not valid:
        return 0.0
    return sum(1 for aa in valid if aa in ACIDIC) / len(valid)


def compute_kr_fraction(seq):
    valid = [aa for aa in seq if aa in STANDARD_AAS]
    if not valid:
        return 0.0
    return sum(1 for aa in valid if aa in BASIC) / len(valid)


def compute_entropy(seq):
    from collections import Counter
    valid = [aa for aa in seq if aa in STANDARD_AAS]
    if len(valid) < 20:
        return None
    counts = Counter(valid)
    total = len(valid)
    ent = 0.0
    for c in counts.values():
        p = c / total
        if p > 0:
            ent -= p * np.log2(p)
    return ent


def compute_disorder_score(seq):
    valid = [aa for aa in seq if aa in STANDARD_AAS]
    if not valid:
        return 0.0
    n_dis = sum(1 for aa in valid if aa in DISORDER_PROMOTING)
    n_ord = sum(1 for aa in valid if aa in ORDER_PROMOTING)
    return (n_dis - n_ord) / len(valid)


def has_de_segment(seq, window=8, threshold=0.375):
    """Check if sequence has a D/E-rich segment."""
    if len(seq) < window:
        return False
    for i in range(len(seq) - window + 1):
        w = seq[i:i+window]
        de = sum(1 for aa in w if aa in ACIDIC)
        if de / window >= threshold:
            return True
    return False


def main():
    base_dir = Path("/storage/kiran-stuff/protein-semantic-collapse")

    # load sequences for flank analysis
    print("Loading full protein sequences...")
    all_seqs = {}
    for fasta in ["sequences/rbp_sequences.fasta", "sequences/nonrbp_sequences.fasta"]:
        for rec in SeqIO.parse(base_dir / fasta, "fasta"):
            parts = rec.id.split("|")
            acc = parts[1] if len(parts) >= 2 else rec.id.split()[0]
            all_seqs[acc] = str(rec.seq)
    print(f"  Loaded {len(all_seqs):,} proteins")

    # load domains
    print("Loading domains...")
    rbd_domains = []
    nonrbd_domains = []

    with open(base_dir / "data/domain_sequences.jsonl") as f:
        for line in f:
            d = json.loads(line)
            if d["is_rbd"]:
                rbd_domains.append(d)
            else:
                nonrbd_domains.append(d)

    print(f"  RBD: {len(rbd_domains):,}, non-RBD: {len(nonrbd_domains):,}")

    # =========================================================
    # ANALYSIS 1: FAMILY OUTLIERS
    # =========================================================
    print("\n" + "=" * 70)
    print("ANALYSIS 1: FAMILY OUTLIERS")
    print("=" * 70)

    # aggregate by family
    family_data = defaultdict(lambda: {
        "de_fracs": [], "kr_fracs": [], "entropies": [],
        "disorder_scores": [], "has_de_segment": [], "flank_de": []
    })

    for d in rbd_domains:
        key = (d["pfam"], d["pfam_name"])
        seq = d["sequence"]

        de = compute_de_fraction(seq)
        kr = compute_kr_fraction(seq)
        ent = compute_entropy(seq)
        dis = compute_disorder_score(seq)
        has_seg = has_de_segment(seq)

        family_data[key]["de_fracs"].append(de)
        family_data[key]["kr_fracs"].append(kr)
        if ent is not None:
            family_data[key]["entropies"].append(ent)
        family_data[key]["disorder_scores"].append(dis)
        family_data[key]["has_de_segment"].append(has_seg)

        # flank D/E
        acc = d["accession"]
        if acc in all_seqs:
            full = all_seqs[acc]
            start, end = d["start"] - 1, d["end"]
            n_flank = full[max(0, start-30):start]
            c_flank = full[end:min(len(full), end+30)]
            flank_de = (compute_de_fraction(n_flank) + compute_de_fraction(c_flank)) / 2
            family_data[key]["flank_de"].append(flank_de)

    # compute family summaries
    family_summary = []
    for (pfam, name), data in family_data.items():
        if len(data["entropies"]) < 10:
            continue
        family_summary.append({
            "pfam": pfam,
            "name": name,
            "n": len(data["de_fracs"]),
            "mean_de": np.mean(data["de_fracs"]),
            "mean_kr": np.mean(data["kr_fracs"]),
            "mean_entropy": np.mean(data["entropies"]),
            "mean_disorder": np.mean(data["disorder_scores"]),
            "pct_with_de_seg": np.mean(data["has_de_segment"]) * 100,
            "mean_flank_de": np.mean(data["flank_de"]) if data["flank_de"] else 0,
        })

    # find outliers: low D/E (unusual for RBDs)
    de_values = [f["mean_de"] for f in family_summary]
    de_mean, de_std = np.mean(de_values), np.std(de_values)

    print(f"\nOverall RBD family D+E: mean={de_mean*100:.2f}%, std={de_std*100:.2f}%")

    print(f"\n{'LOW D+E OUTLIERS (< mean - 1.5*std):'}")
    print(f"{'Family':<12} {'Name':<20} {'N':>5} {'D+E%':>7} {'K+R%':>7} {'Entropy':>8}")
    print("-" * 70)

    low_de_threshold = de_mean - 1.5 * de_std
    low_de_families = [f for f in family_summary if f["mean_de"] < low_de_threshold]
    low_de_families.sort(key=lambda x: x["mean_de"])

    for f in low_de_families[:10]:
        print(f"{f['pfam']:<12} {f['name'][:20]:<20} {f['n']:>5} {f['mean_de']*100:>6.2f}% {f['mean_kr']*100:>6.2f}% {f['mean_entropy']:>8.3f}")

    # high K+R (basic) - potential direct binders without autoinhibition
    kr_values = [f["mean_kr"] for f in family_summary]
    kr_mean, kr_std = np.mean(kr_values), np.std(kr_values)

    print(f"\n{'HIGH K+R OUTLIERS (> mean + 1.5*std) - potential direct binders:'}")
    print(f"{'Family':<12} {'Name':<20} {'N':>5} {'K+R%':>7} {'D+E%':>7} {'Net':>7}")
    print("-" * 70)

    high_kr_threshold = kr_mean + 1.5 * kr_std
    high_kr_families = [f for f in family_summary if f["mean_kr"] > high_kr_threshold]
    high_kr_families.sort(key=lambda x: x["mean_kr"], reverse=True)

    for f in high_kr_families[:10]:
        net = (f["mean_kr"] - f["mean_de"]) * 100
        print(f"{f['pfam']:<12} {f['name'][:20]:<20} {f['n']:>5} {f['mean_kr']*100:>6.2f}% {f['mean_de']*100:>6.2f}% {net:>+6.2f}%")

    # low D/E segment presence
    print(f"\n{'LOW D/E SEGMENT PRESENCE (<30%) - may use different regulation:'}")
    print(f"{'Family':<12} {'Name':<20} {'N':>5} {'%Seg':>7} {'D+E%':>7} {'Disorder':>9}")
    print("-" * 70)

    low_seg_families = [f for f in family_summary if f["pct_with_de_seg"] < 30 and f["n"] >= 20]
    low_seg_families.sort(key=lambda x: x["pct_with_de_seg"])

    for f in low_seg_families[:10]:
        print(f"{f['pfam']:<12} {f['name'][:20]:<20} {f['n']:>5} {f['pct_with_de_seg']:>6.1f}% {f['mean_de']*100:>6.2f}% {f['mean_disorder']:>9.4f}")

    # =========================================================
    # ANALYSIS 2: ENTROPY vs AUTOINHIBITION CORRELATION
    # =========================================================
    print("\n" + "=" * 70)
    print("ANALYSIS 2: ENTROPY vs AUTOINHIBITION CORRELATION")
    print("=" * 70)

    # per-family correlations
    entropies = [f["mean_entropy"] for f in family_summary]
    de_fracs = [f["mean_de"] for f in family_summary]
    de_seg_pcts = [f["pct_with_de_seg"] for f in family_summary]
    flank_des = [f["mean_flank_de"] for f in family_summary if f["mean_flank_de"] > 0]
    entropies_for_flank = [f["mean_entropy"] for f in family_summary if f["mean_flank_de"] > 0]

    # entropy vs domain D/E
    r_de, p_de = stats.pearsonr(entropies, de_fracs)
    print(f"\nEntropy vs Domain D+E:     r={r_de:+.3f}, p={p_de:.2e}")

    # entropy vs D/E segment presence
    r_seg, p_seg = stats.pearsonr(entropies, de_seg_pcts)
    print(f"Entropy vs D/E Seg %:      r={r_seg:+.3f}, p={p_seg:.2e}")

    # entropy vs flank D/E
    if len(flank_des) > 5:
        r_flank, p_flank = stats.pearsonr(entropies_for_flank, flank_des)
        print(f"Entropy vs Flank D+E:      r={r_flank:+.3f}, p={p_flank:.2e}")

    # interpretation
    print(f"\n{'INTERPRETATION:'}")
    if r_de < -0.2 and p_de < 0.05:
        print("  [!] NOVEL: Lower entropy correlates with HIGHER D/E content")
        print("      → Semantic collapse may drive stronger autoinhibition")
    elif r_de > 0.2 and p_de < 0.05:
        print("  [!] Lower entropy correlates with LOWER D/E content")
        print("      → Collapsed vocabulary families use less autoinhibition")
    else:
        print("  [~] No significant correlation between entropy and D/E")
        print("      → Autoinhibition and vocabulary collapse are independent")

    # =========================================================
    # ANALYSIS 3: STRUCTURAL VALIDATION PREP
    # =========================================================
    print("\n" + "=" * 70)
    print("ANALYSIS 3: STRUCTURAL VALIDATION")
    print("=" * 70)

    # identify domains with PDB structures
    # we'll look for accessions that likely have structures
    print("\nIdentifying candidates with likely PDB structures...")

    # for now, output the top outlier families for manual structure lookup
    print(f"\n{'CANDIDATES FOR STRUCTURAL VALIDATION:'}")
    print("(Check if D/E-rich regions show disorder in crystal structures)")
    print()

    # high D/E + high flank D/E families
    high_de_families = sorted(family_summary, key=lambda x: x["mean_de"], reverse=True)[:5]

    print(f"{'High D/E families to validate:'}")
    for f in high_de_families:
        print(f"  {f['pfam']} ({f['name']}): D+E={f['mean_de']*100:.1f}%, n={f['n']}")

    # contrast: low D/E families
    low_de_families_struct = sorted(family_summary, key=lambda x: x["mean_de"])[:5]
    print(f"\n{'Low D/E families (controls):'}")
    for f in low_de_families_struct:
        print(f"  {f['pfam']} ({f['name']}): D+E={f['mean_de']*100:.1f}%, n={f['n']}")

    # =========================================================
    # NOVEL FINDINGS SUMMARY
    # =========================================================
    print("\n" + "=" * 70)
    print("NOVEL FINDINGS SUMMARY")
    print("=" * 70)

    print(f"\n1. OUTLIER FAMILIES:")
    if low_de_families:
        print(f"   - {len([f for f in family_summary if f['mean_de'] < low_de_threshold])} families have unusually LOW D/E")
        print(f"   - Top outlier: {low_de_families[0]['name']} ({low_de_families[0]['mean_de']*100:.1f}% D+E)")

    if high_kr_families:
        print(f"   - {len(high_kr_families)} families have unusually HIGH K+R (direct binders?)")
        print(f"   - Top: {high_kr_families[0]['name']} ({high_kr_families[0]['mean_kr']*100:.1f}% K+R)")

    print(f"\n2. ENTROPY-AUTOINHIBITION LINK:")
    if abs(r_de) > 0.2 and p_de < 0.05:
        direction = "negative" if r_de < 0 else "positive"
        print(f"   - SIGNIFICANT {direction} correlation (r={r_de:.3f})")
    else:
        print(f"   - No significant link (r={r_de:.3f}, p={p_de:.2f})")
        print("   - Semantic collapse and autoinhibition appear INDEPENDENT")


if __name__ == "__main__":
    main()
