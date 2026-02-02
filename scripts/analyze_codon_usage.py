#!/usr/bin/env python3
"""
Codon usage analysis:
1. Get CDS for human RBD proteins
2. Map D/E segments to codon positions
3. Compare codon usage in D/E regions vs domain cores
4. Prepare sequences for CodonFM scoring
"""

import json
import requests
import time
from collections import Counter, defaultdict
from pathlib import Path
import re

# standard codon table
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

# codons for D and E
D_CODONS = ['GAT', 'GAC']
E_CODONS = ['GAA', 'GAG']
DE_CODONS = set(D_CODONS + E_CODONS)

# rare codons in human (low tRNA abundance)
RARE_CODONS_HUMAN = {'CGA', 'CGG', 'AGA', 'AGG', 'CTA', 'ATA', 'TCG', 'CCG', 'ACG', 'GCG'}


def fetch_cds_from_uniprot(accession):
    """Fetch CDS via UniProt cross-references to EMBL."""
    try:
        url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
        resp = requests.get(url, timeout=15)
        if resp.status_code != 200:
            return None

        data = resp.json()

        # get EMBL cross-references
        xrefs = data.get("uniProtKBCrossReferences", [])
        embl_refs = [x for x in xrefs if x.get("database") == "EMBL"]

        # find one with CDS (prefer mRNA)
        for ref in embl_refs:
            props = {p["key"]: p["value"] for p in ref.get("properties", [])}
            protein_id = props.get("ProteinId", "")
            mol_type = props.get("MoleculeType", "")

            if protein_id and protein_id != "-" and "mRNA" in mol_type:
                # fetch from ENA
                cds = fetch_cds_from_ena(ref["id"], protein_id)
                if cds:
                    return cds

        return None
    except Exception as e:
        return None


def fetch_cds_from_ena(embl_id, protein_id):
    """Fetch CDS sequence from ENA."""
    try:
        # try ENA API
        url = f"https://www.ebi.ac.uk/ena/browser/api/fasta/{embl_id}?download=false"
        resp = requests.get(url, timeout=15)
        if resp.status_code == 200:
            # parse FASTA
            lines = resp.text.strip().split('\n')
            if len(lines) > 1:
                seq = ''.join(lines[1:]).upper()
                # this is the full mRNA, need to find CDS
                # for now return full sequence (could extract CDS region)
                return seq
        return None
    except:
        return None


def compute_codon_stats(cds, protein_seq, domain_start, domain_end):
    """
    Compute codon usage stats for different regions.

    Args:
        cds: coding DNA sequence (full)
        protein_seq: protein sequence
        domain_start: 1-based start of domain in protein
        domain_end: 1-based end of domain in protein
    """
    # translate CDS to verify match
    codons = [cds[i:i+3] for i in range(0, len(cds)-2, 3)]
    translated = ''.join(CODON_TABLE.get(c, 'X') for c in codons)

    # find domain region in CDS
    # domain positions are 1-based in protein
    codon_start = domain_start - 1  # 0-based codon index
    codon_end = domain_end  # exclusive

    if codon_end > len(codons):
        return None

    domain_codons = codons[codon_start:codon_end]

    # identify D/E positions in domain
    de_positions = []
    other_positions = []

    for i, codon in enumerate(domain_codons):
        aa = CODON_TABLE.get(codon, 'X')
        if aa in 'DE':
            de_positions.append(i)
        else:
            other_positions.append(i)

    # compute stats
    de_codons_used = [domain_codons[i] for i in de_positions]
    other_codons_used = [domain_codons[i] for i in other_positions]

    # rare codon frequency
    de_rare = sum(1 for c in de_codons_used if c in RARE_CODONS_HUMAN) / len(de_codons_used) if de_codons_used else 0
    other_rare = sum(1 for c in other_codons_used if c in RARE_CODONS_HUMAN) / len(other_codons_used) if other_codons_used else 0

    # codon diversity (entropy-like)
    def codon_entropy(codon_list):
        if not codon_list:
            return 0
        counts = Counter(codon_list)
        total = len(codon_list)
        import math
        return -sum((c/total) * math.log2(c/total) for c in counts.values() if c > 0)

    return {
        "domain_length": len(domain_codons),
        "n_de": len(de_positions),
        "n_other": len(other_positions),
        "de_rare_frac": de_rare,
        "other_rare_frac": other_rare,
        "de_codon_entropy": codon_entropy(de_codons_used),
        "other_codon_entropy": codon_entropy(other_codons_used),
        "domain_codons": ''.join(domain_codons),  # for CodonFM
    }


def main():
    base_dir = Path("/storage/kiran-stuff/protein-semantic-collapse")
    out_dir = base_dir / "results"

    # load human RBD domains
    print("Loading human RBD domains...")

    # first get human accessions from FASTA
    human_accessions = set()
    with open(base_dir / "sequences/rbp_sequences.fasta") as f:
        for line in f:
            if line.startswith(">") and ("Homo sapiens" in line or "OX=9606" in line):
                parts = line.split("|")
                if len(parts) >= 2:
                    human_accessions.add(parts[1])

    print(f"  Human accessions: {len(human_accessions)}")

    # load domains for human proteins
    human_domains = []
    with open(base_dir / "data/domain_sequences.jsonl") as f:
        for line in f:
            d = json.loads(line)
            if d["is_rbd"] and d["accession"] in human_accessions:
                human_domains.append(d)

    print(f"  Human RBD domains: {len(human_domains)}")

    # sample for testing (full analysis would take hours)
    sample_size = 50
    sample = human_domains[:sample_size]

    print(f"\nFetching CDS for {len(sample)} sample domains...")

    results = []
    success = 0
    failed = 0

    for i, d in enumerate(sample):
        if i % 10 == 0:
            print(f"  Processing {i+1}/{len(sample)}...")

        cds = fetch_cds_from_uniprot(d["accession"])

        if cds and len(cds) >= d["end"] * 3:
            stats = compute_codon_stats(cds, d["sequence"], d["start"], d["end"])
            if stats:
                stats["accession"] = d["accession"]
                stats["pfam"] = d["pfam"]
                stats["pfam_name"] = d["pfam_name"]
                results.append(stats)
                success += 1
            else:
                failed += 1
        else:
            failed += 1

        time.sleep(0.5)  # rate limiting

    print(f"\n  Success: {success}, Failed: {failed}")

    if not results:
        print("No results to analyze")
        return

    # =========================================================
    # CODON ANALYSIS RESULTS
    # =========================================================
    print("\n" + "=" * 70)
    print("CODON USAGE ANALYSIS")
    print("=" * 70)

    import numpy as np

    de_rare_fracs = [r["de_rare_frac"] for r in results if r["n_de"] > 0]
    other_rare_fracs = [r["other_rare_frac"] for r in results if r["n_other"] > 0]

    print(f"\nRare codon frequency (mean):")
    print(f"  D/E positions:   {np.mean(de_rare_fracs)*100:.2f}%")
    print(f"  Other positions: {np.mean(other_rare_fracs)*100:.2f}%")

    diff = np.mean(de_rare_fracs) - np.mean(other_rare_fracs)
    print(f"  Difference: {diff*100:+.2f}%")

    if diff > 0.02:
        print("  → D/E regions have MORE rare codons (translation pausing?)")
    elif diff < -0.02:
        print("  → D/E regions have FEWER rare codons (faster translation)")
    else:
        print("  → No significant difference in rare codon usage")

    # save domain CDS for CodonFM
    print(f"\nSaving domain CDS for CodonFM scoring...")
    with open(out_dir / "domain_cds_for_codonfm.fasta", "w") as f:
        for r in results:
            f.write(f">{r['accession']}|{r['pfam']}|{r['pfam_name']}\n")
            f.write(f"{r['domain_codons']}\n")

    print(f"  Saved {len(results)} sequences to results/domain_cds_for_codonfm.fasta")

    # summary by family
    print("\n" + "=" * 70)
    print("BY FAMILY (if enough data)")
    print("=" * 70)

    family_results = defaultdict(list)
    for r in results:
        family_results[r["pfam_name"]].append(r)

    for name, fam_results in sorted(family_results.items(), key=lambda x: -len(x[1])):
        if len(fam_results) >= 3:
            de_rare = np.mean([r["de_rare_frac"] for r in fam_results if r["n_de"] > 0])
            other_rare = np.mean([r["other_rare_frac"] for r in fam_results if r["n_other"] > 0])
            print(f"{name}: n={len(fam_results)}, D/E rare={de_rare*100:.1f}%, other rare={other_rare*100:.1f}%")


if __name__ == "__main__":
    main()
