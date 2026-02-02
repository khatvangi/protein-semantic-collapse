#!/usr/bin/env python3
"""
AlphaFold pLDDT validation:
- Fetch AlphaFold structures for RBD proteins
- Extract pLDDT scores (confidence = order proxy)
- Compare pLDDT in D/E regions vs domain cores
"""

import json
import requests
import numpy as np
from collections import defaultdict
from pathlib import Path
import time

PLDDT_DISORDERED_THRESHOLD = 70


def fetch_alphafold_plddt(accession):
    """Fetch per-residue pLDDT scores from AlphaFold DB."""
    try:
        # get prediction metadata
        url = f"https://alphafold.ebi.ac.uk/api/prediction/{accession}"
        resp = requests.get(url, timeout=30)
        if resp.status_code != 200:
            return None

        data = resp.json()
        if not data:
            return None

        entry = data[0] if isinstance(data, list) else data

        # get pLDDT document URL
        plddt_url = entry.get('plddtDocUrl')
        if not plddt_url:
            return None

        # fetch per-residue pLDDT
        resp2 = requests.get(plddt_url, timeout=30)
        if resp2.status_code != 200:
            return None

        plddt_data = resp2.json()

        # extract scores
        if isinstance(plddt_data, dict):
            scores = plddt_data.get('confidenceScore', [])
            return scores if scores else None

        return None

    except Exception as e:
        return None


def analyze_plddt_by_region(plddt_scores, sequence, domain_start, domain_end):
    """Analyze pLDDT scores for different regions."""
    if not plddt_scores or len(plddt_scores) < domain_end:
        return None

    # domain region (0-indexed)
    domain_start_idx = domain_start - 1
    domain_end_idx = domain_end

    domain_plddt = plddt_scores[domain_start_idx:domain_end_idx]
    domain_seq = sequence

    if len(domain_plddt) != len(domain_seq):
        min_len = min(len(domain_plddt), len(domain_seq))
        domain_plddt = domain_plddt[:min_len]
        domain_seq = domain_seq[:min_len]

    # separate D/E vs other
    de_plddt = []
    other_plddt = []

    for aa, score in zip(domain_seq, domain_plddt):
        if aa in 'DE':
            de_plddt.append(score)
        else:
            other_plddt.append(score)

    # flank pLDDT
    flank_size = 30
    n_flank_start = max(0, domain_start_idx - flank_size)
    n_flank_plddt = plddt_scores[n_flank_start:domain_start_idx]

    c_flank_end = min(len(plddt_scores), domain_end_idx + flank_size)
    c_flank_plddt = plddt_scores[domain_end_idx:c_flank_end]

    return {
        "mean_de": np.mean(de_plddt) if de_plddt else None,
        "mean_other": np.mean(other_plddt) if other_plddt else None,
        "mean_domain": np.mean(domain_plddt) if domain_plddt else None,
        "mean_n_flank": np.mean(n_flank_plddt) if n_flank_plddt else None,
        "mean_c_flank": np.mean(c_flank_plddt) if c_flank_plddt else None,
        "n_de": len(de_plddt),
        "n_other": len(other_plddt),
        "pct_de_disordered": sum(1 for p in de_plddt if p < PLDDT_DISORDERED_THRESHOLD) / len(de_plddt) * 100 if de_plddt else 0,
        "pct_flank_disordered": sum(1 for p in (n_flank_plddt + c_flank_plddt) if p < PLDDT_DISORDERED_THRESHOLD) / len(n_flank_plddt + c_flank_plddt) * 100 if (n_flank_plddt + c_flank_plddt) else 0,
    }


def main():
    base_dir = Path("/storage/kiran-stuff/protein-semantic-collapse")

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

    print(f"  RBD: {len(rbd_domains)}, non-RBD: {len(nonrbd_domains)}")

    # sample
    sample_size = 100
    rbd_sample = rbd_domains[:sample_size]
    nonrbd_sample = nonrbd_domains[:sample_size]

    print(f"\nFetching AlphaFold pLDDT...")

    rbd_results = []
    nonrbd_results = []

    # RBD
    print("\nRBD domains:")
    for i, d in enumerate(rbd_sample):
        if i % 20 == 0:
            print(f"  {i+1}/{len(rbd_sample)}... ({len(rbd_results)} success)")

        plddt = fetch_alphafold_plddt(d["accession"])
        if plddt:
            analysis = analyze_plddt_by_region(plddt, d["sequence"], d["start"], d["end"])
            if analysis and analysis["mean_de"] is not None:
                analysis["accession"] = d["accession"]
                analysis["pfam"] = d["pfam"]
                rbd_results.append(analysis)

        time.sleep(0.2)

    # non-RBD
    print("\nnon-RBD domains:")
    for i, d in enumerate(nonrbd_sample):
        if i % 20 == 0:
            print(f"  {i+1}/{len(nonrbd_sample)}... ({len(nonrbd_results)} success)")

        plddt = fetch_alphafold_plddt(d["accession"])
        if plddt:
            analysis = analyze_plddt_by_region(plddt, d["sequence"], d["start"], d["end"])
            if analysis and analysis["mean_de"] is not None:
                analysis["accession"] = d["accession"]
                analysis["pfam"] = d["pfam"]
                nonrbd_results.append(analysis)

        time.sleep(0.2)

    print(f"\n  RBD with pLDDT: {len(rbd_results)}")
    print(f"  non-RBD with pLDDT: {len(nonrbd_results)}")

    if not rbd_results:
        print("No results")
        return

    # =========================================================
    # RESULTS
    # =========================================================
    print("\n" + "=" * 70)
    print("ALPHAFOLD pLDDT ANALYSIS")
    print("=" * 70)

    def summarize(results, label):
        de = [r["mean_de"] for r in results if r["mean_de"]]
        other = [r["mean_other"] for r in results if r["mean_other"]]
        domain = [r["mean_domain"] for r in results if r["mean_domain"]]
        nflank = [r["mean_n_flank"] for r in results if r["mean_n_flank"]]
        cflank = [r["mean_c_flank"] for r in results if r["mean_c_flank"]]
        de_dis = [r["pct_de_disordered"] for r in results]
        flank_dis = [r["pct_flank_disordered"] for r in results]

        print(f"\n{label} (n={len(results)}):")
        print(f"  {'Region':<25} {'Mean pLDDT':>12} {'% Disordered':>15}")
        print(f"  {'-'*55}")
        print(f"  {'D/E positions':<25} {np.mean(de):>12.1f} {np.mean(de_dis):>14.1f}%")
        print(f"  {'Other positions':<25} {np.mean(other):>12.1f}")
        print(f"  {'Domain overall':<25} {np.mean(domain):>12.1f}")
        print(f"  {'N-terminal flank':<25} {np.mean(nflank):>12.1f} {np.mean(flank_dis):>14.1f}%")
        print(f"  {'C-terminal flank':<25} {np.mean(cflank):>12.1f}")

        return {
            "de": np.mean(de), "other": np.mean(other), "domain": np.mean(domain),
            "nflank": np.mean(nflank), "cflank": np.mean(cflank),
            "de_dis": np.mean(de_dis), "flank_dis": np.mean(flank_dis)
        }

    rbd_stats = summarize(rbd_results, "RBD DOMAINS")
    nonrbd_stats = summarize(nonrbd_results, "non-RBD DOMAINS") if nonrbd_results else None

    # =========================================================
    # COMPARISON
    # =========================================================
    print("\n" + "=" * 70)
    print("KEY COMPARISONS")
    print("=" * 70)

    de_vs_other = rbd_stats["de"] - rbd_stats["other"]
    flank_vs_domain = (rbd_stats["nflank"] + rbd_stats["cflank"]) / 2 - rbd_stats["domain"]

    print(f"\n1. D/E vs other positions (RBD): {de_vs_other:+.1f} pLDDT")
    if de_vs_other < -3:
        print("   ✓ D/E positions MORE disordered - SUPPORTS hypothesis")
    elif de_vs_other > 3:
        print("   ✗ D/E positions MORE ordered - unexpected")
    else:
        print("   ~ No significant difference")

    print(f"\n2. Flanks vs domain core (RBD): {flank_vs_domain:+.1f} pLDDT")
    if flank_vs_domain < -5:
        print("   ✓ Flanks MORE disordered - SUPPORTS IDR hypothesis")
    elif flank_vs_domain > 5:
        print("   ✗ Flanks MORE ordered - unexpected")
    else:
        print("   ~ Similar disorder levels")

    if nonrbd_stats:
        rbd_vs_nonrbd_de = rbd_stats["de"] - nonrbd_stats["de"]
        print(f"\n3. RBD vs non-RBD D/E pLDDT: {rbd_vs_nonrbd_de:+.1f}")
        if rbd_vs_nonrbd_de < -3:
            print("   ✓ RBD D/E positions MORE disordered than non-RBD")
        elif rbd_vs_nonrbd_de > 3:
            print("   ✗ RBD D/E positions MORE ordered")
        else:
            print("   ~ Similar")

    print(f"\n4. % disordered (pLDDT<70):")
    print(f"   RBD D/E positions: {rbd_stats['de_dis']:.1f}%")
    print(f"   RBD flanks: {rbd_stats['flank_dis']:.1f}%")


if __name__ == "__main__":
    main()
