#!/usr/bin/env python3
"""
Functional analysis:
1. GO enrichment for RBD families (high D/E vs low D/E outliers)
2. Disease association check via UniProt
"""

import json
import gzip
from collections import defaultdict, Counter
from pathlib import Path
import requests

# GO term categories of interest
GO_CATEGORIES = {
    "GO:0003723": "RNA binding",
    "GO:0003676": "nucleic acid binding",
    "GO:0005840": "ribosome",
    "GO:0006412": "translation",
    "GO:0000398": "mRNA splicing",
    "GO:0006397": "mRNA processing",
    "GO:0016070": "RNA metabolic process",
    "GO:0010467": "gene expression",
    "GO:0005634": "nucleus",
    "GO:0005737": "cytoplasm",
    "GO:0005829": "cytosol",
    "GO:0005739": "mitochondrion",
    "GO:0006950": "response to stress",
}


def load_go_annotations(go_file, accessions):
    """Load GO annotations for specified accessions."""
    acc_go = defaultdict(set)

    with gzip.open(go_file, 'rt') as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split(',')
            if len(parts) >= 6:
                acc = parts[2]  # SP_PRIMARY
                go_id = parts[5]  # GO_ID
                if acc in accessions:
                    acc_go[acc].add(go_id)

    return acc_go


def get_disease_annotations(accession):
    """Fetch disease annotations from UniProt for an accession."""
    try:
        url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
        resp = requests.get(url, timeout=10)
        if resp.status_code != 200:
            return []

        data = resp.json()
        diseases = []

        # check disease involvement
        for comment in data.get("comments", []):
            if comment.get("commentType") == "DISEASE":
                disease = comment.get("disease", {})
                if disease:
                    diseases.append({
                        "name": disease.get("diseaseId", ""),
                        "description": disease.get("description", ""),
                        "acronym": disease.get("acronym", "")
                    })

        return diseases
    except Exception as e:
        return []


def main():
    base_dir = Path("/storage/kiran-stuff/protein-semantic-collapse")

    # load domain data
    print("Loading domain data...")
    family_proteins = defaultdict(set)  # pfam -> set of accessions
    family_info = {}  # pfam -> (name, is_rbd)

    # track high D/E and low D/E families
    high_de_families = ["PF00658", "PF00012", "PF00076"]  # PABP, HSP70, RRM_1
    low_de_families = ["PF00189", "PF01479", "PF01280"]   # Ribosomal

    with open(base_dir / "data/domain_sequences.jsonl") as f:
        for line in f:
            d = json.loads(line)
            if d["is_rbd"]:
                family_proteins[d["pfam"]].add(d["accession"])
                family_info[d["pfam"]] = (d["pfam_name"], d["is_rbd"])

    # get all accessions
    all_rbd_accessions = set()
    for accs in family_proteins.values():
        all_rbd_accessions.update(accs)

    print(f"  Total RBD families: {len(family_proteins)}")
    print(f"  Total unique accessions: {len(all_rbd_accessions)}")

    # load GO annotations
    print("\nLoading GO annotations...")
    go_file = base_dir / "data/clean_pdb_chain_go.csv.gz"
    acc_go = load_go_annotations(go_file, all_rbd_accessions)
    print(f"  Proteins with GO: {len(acc_go)}")

    # =========================================================
    # GO ENRICHMENT COMPARISON
    # =========================================================
    print("\n" + "=" * 70)
    print("GO TERM COMPARISON: HIGH D/E vs LOW D/E FAMILIES")
    print("=" * 70)

    # collect GO terms by group
    high_de_go = Counter()
    low_de_go = Counter()
    high_de_n = 0
    low_de_n = 0

    for pfam in high_de_families:
        for acc in family_proteins.get(pfam, []):
            if acc in acc_go:
                high_de_n += 1
                high_de_go.update(acc_go[acc])

    for pfam in low_de_families:
        for acc in family_proteins.get(pfam, []):
            if acc in acc_go:
                low_de_n += 1
                low_de_go.update(acc_go[acc])

    print(f"\nHigh D/E group: {high_de_n} proteins with GO annotations")
    print(f"Low D/E group: {low_de_n} proteins with GO annotations")

    # compare key GO terms
    print(f"\n{'GO Term':<50} {'High D/E':>10} {'Low D/E':>10} {'Diff':>10}")
    print("-" * 80)

    for go_id, name in GO_CATEGORIES.items():
        high_count = high_de_go.get(go_id, 0)
        low_count = low_de_go.get(go_id, 0)

        high_pct = (high_count / high_de_n * 100) if high_de_n > 0 else 0
        low_pct = (low_count / low_de_n * 100) if low_de_n > 0 else 0
        diff = high_pct - low_pct

        if high_count > 0 or low_count > 0:
            marker = "***" if abs(diff) > 20 else "**" if abs(diff) > 10 else ""
            print(f"{name[:48]:<50} {high_pct:>9.1f}% {low_pct:>9.1f}% {diff:>+9.1f}% {marker}")

    # top enriched in each group
    print(f"\n{'TOP GO TERMS UNIQUE TO HIGH D/E:'}")
    high_only = [(go, c) for go, c in high_de_go.most_common(20) if low_de_go.get(go, 0) == 0]
    for go, count in high_only[:5]:
        pct = count / high_de_n * 100
        print(f"  {go}: {count} ({pct:.1f}%)")

    print(f"\n{'TOP GO TERMS UNIQUE TO LOW D/E (ribosomal):'}")
    low_only = [(go, c) for go, c in low_de_go.most_common(20) if high_de_go.get(go, 0) == 0]
    for go, count in low_only[:5]:
        pct = count / low_de_n * 100
        print(f"  {go}: {count} ({pct:.1f}%)")

    # =========================================================
    # DISEASE ASSOCIATION CHECK
    # =========================================================
    print("\n" + "=" * 70)
    print("DISEASE ASSOCIATIONS (sample check)")
    print("=" * 70)

    # sample a few proteins from each group
    print("\nChecking disease associations for sample proteins...")

    high_de_diseases = []
    low_de_diseases = []

    # sample from high D/E
    sample_high = []
    for pfam in high_de_families:
        sample_high.extend(list(family_proteins.get(pfam, []))[:3])

    print(f"\nHigh D/E family samples ({len(sample_high[:9])} proteins):")
    for acc in sample_high[:9]:
        diseases = get_disease_annotations(acc)
        if diseases:
            high_de_diseases.extend(diseases)
            print(f"  {acc}: {len(diseases)} disease(s) - {diseases[0]['name'][:40] if diseases else ''}")

    # sample from low D/E
    sample_low = []
    for pfam in low_de_families:
        sample_low.extend(list(family_proteins.get(pfam, []))[:3])

    print(f"\nLow D/E family samples ({len(sample_low[:9])} proteins):")
    for acc in sample_low[:9]:
        diseases = get_disease_annotations(acc)
        if diseases:
            low_de_diseases.extend(diseases)
            print(f"  {acc}: {len(diseases)} disease(s) - {diseases[0]['name'][:40] if diseases else ''}")

    # summary
    print("\n" + "=" * 70)
    print("FUNCTIONAL ANALYSIS SUMMARY")
    print("=" * 70)

    print(f"\n1. GO ENRICHMENT:")
    # find biggest differences
    all_go = set(high_de_go.keys()) | set(low_de_go.keys())
    diffs = []
    for go in all_go:
        h = (high_de_go.get(go, 0) / high_de_n * 100) if high_de_n > 0 else 0
        l = (low_de_go.get(go, 0) / low_de_n * 100) if low_de_n > 0 else 0
        if h > 5 or l > 5:  # only consider reasonably common terms
            diffs.append((go, h - l, h, l))

    diffs.sort(key=lambda x: x[1], reverse=True)

    print("   Enriched in HIGH D/E (non-ribosomal RBDs):")
    for go, diff, h, l in diffs[:3]:
        name = GO_CATEGORIES.get(go, go)
        print(f"     {name}: {h:.1f}% vs {l:.1f}% ({diff:+.1f}%)")

    print("   Enriched in LOW D/E (ribosomal proteins):")
    for go, diff, h, l in diffs[-3:]:
        name = GO_CATEGORIES.get(go, go)
        print(f"     {name}: {l:.1f}% vs {h:.1f}% ({-diff:+.1f}%)")

    print(f"\n2. DISEASE ASSOCIATIONS:")
    print(f"   High D/E samples: {len(high_de_diseases)} diseases found")
    print(f"   Low D/E samples: {len(low_de_diseases)} diseases found")


if __name__ == "__main__":
    main()
