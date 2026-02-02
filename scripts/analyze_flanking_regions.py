#!/usr/bin/env python3
"""
Analyze flanking regions around domains to test if D/E-rich IDRs
lie outside Pfam domain boundaries.
"""

import json
import numpy as np
from collections import Counter
from pathlib import Path
from Bio import SeqIO

STANDARD_AAS = set("ACDEFGHIKLMNPQRSTVWY")
ACIDIC = set("DE")
DISORDER_PROMOTING = set("AEGRQSPK")
ORDER_PROMOTING = set("WFYILMVNC")

FLANK_SIZE = 30  # residues on each side


def load_sequences(fasta_file):
    """Load sequences into dict {accession: sequence}."""
    seqs = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        # extract accession from header (sp|ACCESSION|NAME format)
        parts = record.id.split("|")
        if len(parts) >= 2:
            acc = parts[1]
        else:
            acc = record.id.split()[0]
        seqs[acc] = str(record.seq)
    return seqs


def compute_de_fraction(seq):
    """Fraction of D+E in sequence."""
    if not seq:
        return 0.0
    valid = [aa for aa in seq if aa in STANDARD_AAS]
    if not valid:
        return 0.0
    de_count = sum(1 for aa in valid if aa in ACIDIC)
    return de_count / len(valid)


def compute_disorder_score(seq):
    """Simple disorder score based on composition."""
    if not seq or len(seq) < 5:
        return 0.0
    valid = [aa for aa in seq if aa in STANDARD_AAS]
    if not valid:
        return 0.0
    n_dis = sum(1 for aa in valid if aa in DISORDER_PROMOTING)
    n_ord = sum(1 for aa in valid if aa in ORDER_PROMOTING)
    return (n_dis - n_ord) / len(valid)


def analyze_flanks(domains, full_seqs, label):
    """Analyze flanking regions for a set of domains."""
    results = {
        "n_total": len(domains),
        "n_with_flanks": 0,
        "n_term_de": [],
        "c_term_de": [],
        "domain_de": [],
        "n_term_disorder": [],
        "c_term_disorder": [],
        "domain_disorder": [],
    }

    for d in domains:
        acc = d["accession"]
        if acc not in full_seqs:
            continue

        full_seq = full_seqs[acc]
        start = d["start"] - 1  # convert to 0-indexed
        end = d["end"]

        # extract flanks
        n_flank_start = max(0, start - FLANK_SIZE)
        n_flank = full_seq[n_flank_start:start]

        c_flank_end = min(len(full_seq), end + FLANK_SIZE)
        c_flank = full_seq[end:c_flank_end]

        domain_seq = d["sequence"]

        # only count if we have meaningful flanks
        if len(n_flank) >= 10 or len(c_flank) >= 10:
            results["n_with_flanks"] += 1

            if len(n_flank) >= 10:
                results["n_term_de"].append(compute_de_fraction(n_flank))
                results["n_term_disorder"].append(compute_disorder_score(n_flank))

            if len(c_flank) >= 10:
                results["c_term_de"].append(compute_de_fraction(c_flank))
                results["c_term_disorder"].append(compute_disorder_score(c_flank))

            results["domain_de"].append(compute_de_fraction(domain_seq))
            results["domain_disorder"].append(compute_disorder_score(domain_seq))

    return results


def main():
    base_dir = Path("/storage/kiran-stuff/protein-semantic-collapse")

    # load full sequences
    print("Loading full protein sequences...")
    rbp_seqs = load_sequences(base_dir / "sequences/rbp_sequences.fasta")
    nonrbp_seqs = load_sequences(base_dir / "sequences/nonrbp_sequences.fasta")
    all_seqs = {**rbp_seqs, **nonrbp_seqs}
    print(f"  Loaded {len(all_seqs):,} proteins")

    # load domains
    print("Loading domain data...")
    rbd_domains = []
    nonrbd_domains = []

    with open(base_dir / "data/domain_sequences.jsonl") as f:
        for line in f:
            d = json.loads(line)
            if d["is_rbd"]:
                rbd_domains.append(d)
            else:
                nonrbd_domains.append(d)

    print(f"  RBD domains: {len(rbd_domains):,}")
    print(f"  non-RBD domains: {len(nonrbd_domains):,}")

    # analyze flanks
    print(f"\nAnalyzing flanking regions (±{FLANK_SIZE} residues)...")
    rbd_flanks = analyze_flanks(rbd_domains, all_seqs, "RBD")
    nonrbd_flanks = analyze_flanks(nonrbd_domains, all_seqs, "non-RBD")

    # report
    print("\n" + "=" * 70)
    print("FLANKING REGION ANALYSIS")
    print("=" * 70)

    print(f"\n{'Metric':<40} {'RBD':>12} {'non-RBD':>12}")
    print("-" * 65)
    print(f"{'Domains with flanks':<40} {rbd_flanks['n_with_flanks']:>12,} {nonrbd_flanks['n_with_flanks']:>12,}")

    # D/E content
    print(f"\n{'D+E Fraction (mean):'}")
    print(f"  {'N-terminal flank':<36} {np.mean(rbd_flanks['n_term_de'])*100:>11.2f}% {np.mean(nonrbd_flanks['n_term_de'])*100:>11.2f}%")
    print(f"  {'Domain (for comparison)':<36} {np.mean(rbd_flanks['domain_de'])*100:>11.2f}% {np.mean(nonrbd_flanks['domain_de'])*100:>11.2f}%")
    print(f"  {'C-terminal flank':<36} {np.mean(rbd_flanks['c_term_de'])*100:>11.2f}% {np.mean(nonrbd_flanks['c_term_de'])*100:>11.2f}%")

    # highlight differences
    rbd_n_de = np.mean(rbd_flanks['n_term_de']) * 100
    rbd_c_de = np.mean(rbd_flanks['c_term_de']) * 100
    rbd_dom_de = np.mean(rbd_flanks['domain_de']) * 100

    nonrbd_n_de = np.mean(nonrbd_flanks['n_term_de']) * 100
    nonrbd_c_de = np.mean(nonrbd_flanks['c_term_de']) * 100

    print(f"\n  RBD N-flank vs domain: {rbd_n_de - rbd_dom_de:+.2f}%")
    print(f"  RBD C-flank vs domain: {rbd_c_de - rbd_dom_de:+.2f}%")

    # disorder scores
    print(f"\n{'Disorder Score (mean):'}")
    print(f"  {'N-terminal flank':<36} {np.mean(rbd_flanks['n_term_disorder']):>12.4f} {np.mean(nonrbd_flanks['n_term_disorder']):>12.4f}")
    print(f"  {'Domain':<36} {np.mean(rbd_flanks['domain_disorder']):>12.4f} {np.mean(nonrbd_flanks['domain_disorder']):>12.4f}")
    print(f"  {'C-terminal flank':<36} {np.mean(rbd_flanks['c_term_disorder']):>12.4f} {np.mean(nonrbd_flanks['c_term_disorder']):>12.4f}")

    # key comparison
    rbd_n_dis = np.mean(rbd_flanks['n_term_disorder'])
    rbd_c_dis = np.mean(rbd_flanks['c_term_disorder'])
    rbd_dom_dis = np.mean(rbd_flanks['domain_disorder'])

    print(f"\n  RBD N-flank vs domain: {rbd_n_dis - rbd_dom_dis:+.4f}")
    print(f"  RBD C-flank vs domain: {rbd_c_dis - rbd_dom_dis:+.4f}")

    # summary
    print("\n" + "=" * 70)
    print("INTERPRETATION")
    print("=" * 70)

    # test: are flanks more disordered than domains?
    flank_more_disordered = (rbd_n_dis > rbd_dom_dis) or (rbd_c_dis > rbd_dom_dis)

    # test: are flanks more acidic than domains?
    flank_more_acidic = (rbd_n_de > rbd_dom_de) or (rbd_c_de > rbd_dom_de)

    if flank_more_disordered:
        print("\n[+] Flanking regions ARE more disordered than domain cores")
        print("    → Supports hypothesis that IDRs lie outside Pfam boundaries")
    else:
        print("\n[-] Flanking regions are NOT more disordered")

    if flank_more_acidic:
        which = []
        if rbd_n_de > rbd_dom_de:
            which.append("N-terminal")
        if rbd_c_de > rbd_dom_de:
            which.append("C-terminal")
        print(f"\n[+] {' and '.join(which)} flanks ARE more acidic than domains")
        print("    → D/E-rich autoinhibitory regions may extend beyond domain")
    else:
        print("\n[-] Flanks are NOT more acidic than domains")

    # RBD vs non-RBD flank comparison
    print(f"\n{'RBD vs non-RBD flank comparison:'}")
    n_diff = rbd_n_de - nonrbd_n_de
    c_diff = rbd_c_de - nonrbd_c_de
    print(f"  N-flank D+E: RBD {rbd_n_de:.2f}% vs non-RBD {nonrbd_n_de:.2f}% ({n_diff:+.2f}%)")
    print(f"  C-flank D+E: RBD {rbd_c_de:.2f}% vs non-RBD {nonrbd_c_de:.2f}% ({c_diff:+.2f}%)")


if __name__ == "__main__":
    main()
