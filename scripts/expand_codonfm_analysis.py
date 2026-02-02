#!/usr/bin/env python3
"""
Expand CodonFM analysis to 100+ properly curated RBD domains.

The issue with the previous analysis was that CDS extraction had errors
(internal stop codons). This script:
1. Fetches CDS from UniProt/Ensembl
2. Validates coding sequence integrity
3. Extracts domain-specific coding sequences
4. Runs CodonFM scoring
"""

import json
import requests
import time
from pathlib import Path
from collections import defaultdict
import subprocess
import os

BASE_DIR = Path("/storage/kiran-stuff/protein-semantic-collapse")
DATA_DIR = BASE_DIR / "data"
RESULTS_DIR = BASE_DIR / "results"
CODONFM_DIR = RESULTS_DIR / "codonfm_analysis"

# codon tables
CODON_TO_AA = {
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

# D/E codons
DE_CODONS = {'GAT', 'GAC', 'GAA', 'GAG'}


def translate(cds):
    """translate CDS to protein"""
    protein = []
    for i in range(0, len(cds) - 2, 3):
        codon = cds[i:i+3].upper()
        aa = CODON_TO_AA.get(codon, 'X')
        protein.append(aa)
    return ''.join(protein)


def validate_cds(cds, expected_protein):
    """check if CDS translates to expected protein"""
    # remove stop codon if present
    if len(cds) % 3 != 0:
        return False, "Length not divisible by 3"

    translated = translate(cds)

    # check for internal stops
    if '*' in translated[:-1]:  # internal stop
        return False, f"Internal stop codon at position {translated.find('*')}"

    # check match (allowing for some flexibility)
    if translated.rstrip('*') == expected_protein:
        return True, "Perfect match"

    # check if close enough (>95% identity)
    matches = sum(1 for a, b in zip(translated, expected_protein) if a == b)
    identity = matches / max(len(translated), len(expected_protein))

    if identity > 0.95:
        return True, f"Match with {identity*100:.1f}% identity"

    return False, f"Poor match ({identity*100:.1f}% identity)"


def fetch_cds_from_uniprot(accession):
    """fetch CDS sequence from UniProt cross-references"""
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"

    try:
        resp = requests.get(url, timeout=30)
        if resp.status_code != 200:
            return None

        data = resp.json()

        # get protein sequence for validation
        protein_seq = data.get('sequence', {}).get('value', '')

        # look for EMBL/RefSeq cross-references with CDS
        cds_refs = []
        for xref in data.get('uniProtKBCrossReferences', []):
            if xref.get('database') in ['EMBL', 'RefSeq']:
                for prop in xref.get('properties', []):
                    if prop.get('key') == 'ProteinId':
                        cds_refs.append({
                            'database': xref['database'],
                            'nucleotide_id': xref.get('id'),
                            'protein_id': prop.get('value')
                        })

        return {
            'accession': accession,
            'protein_sequence': protein_seq,
            'cds_refs': cds_refs
        }

    except Exception as e:
        print(f"Error fetching {accession}: {e}")
        return None


def fetch_cds_sequence(embl_id):
    """fetch CDS sequence from ENA/EMBL"""
    url = f"https://www.ebi.ac.uk/ena/browser/api/fasta/{embl_id}?download=true&gzip=false"

    try:
        resp = requests.get(url, timeout=30)
        if resp.status_code == 200:
            lines = resp.text.strip().split('\n')
            seq = ''.join(line for line in lines if not line.startswith('>'))
            return seq.upper().replace('U', 'T')
        return None
    except:
        return None


def extract_domain_cds(full_cds, protein_start, protein_end):
    """extract CDS for domain region (1-indexed protein positions)"""
    # convert protein positions to CDS positions
    cds_start = (protein_start - 1) * 3  # 0-indexed
    cds_end = protein_end * 3

    if cds_end > len(full_cds):
        return None

    return full_cds[cds_start:cds_end]


def compute_codon_nll_simple(cds, model_scores=None):
    """
    Simple codon analysis without CodonFM model.

    Uses codon adaptation index (CAI) as proxy for optimization.
    Human codon usage table from Kazusa database.
    """
    # human codon usage frequencies (per 1000)
    HUMAN_CODON_FREQ = {
        'TTT': 17.6, 'TTC': 20.3, 'TTA': 7.7, 'TTG': 12.9,
        'TCT': 15.2, 'TCC': 17.7, 'TCA': 12.2, 'TCG': 4.4,
        'TAT': 12.2, 'TAC': 15.3, 'TAA': 1.0, 'TAG': 0.8,
        'TGT': 10.6, 'TGC': 12.6, 'TGA': 1.6, 'TGG': 13.2,
        'CTT': 13.2, 'CTC': 19.6, 'CTA': 7.2, 'CTG': 39.6,
        'CCT': 17.5, 'CCC': 19.8, 'CCA': 16.9, 'CCG': 6.9,
        'CAT': 10.9, 'CAC': 15.1, 'CAA': 12.3, 'CAG': 34.2,
        'CGT': 4.5, 'CGC': 10.4, 'CGA': 6.2, 'CGG': 11.4,
        'ATT': 16.0, 'ATC': 20.8, 'ATA': 7.5, 'ATG': 22.0,
        'ACT': 13.1, 'ACC': 18.9, 'ACA': 15.1, 'ACG': 6.1,
        'AAT': 17.0, 'AAC': 19.1, 'AAA': 24.4, 'AAG': 31.9,
        'AGT': 12.1, 'AGC': 19.5, 'AGA': 12.2, 'AGG': 12.0,
        'GTT': 11.0, 'GTC': 14.5, 'GTA': 7.1, 'GTG': 28.1,
        'GCT': 18.4, 'GCC': 27.7, 'GCA': 15.8, 'GCG': 7.4,
        'GAT': 21.8, 'GAC': 25.1, 'GAA': 29.0, 'GAG': 39.6,
        'GGT': 10.8, 'GGC': 22.2, 'GGA': 16.5, 'GGG': 16.5,
    }

    # compute per-position CAI-like score (negative log frequency)
    scores = []
    for i in range(0, len(cds) - 2, 3):
        codon = cds[i:i+3].upper()
        if codon in HUMAN_CODON_FREQ:
            freq = HUMAN_CODON_FREQ[codon]
            # negative log frequency (higher = rarer)
            score = -1 * (freq / 1000)  # normalize
            scores.append((i // 3, codon, score))

    return scores


def analyze_de_vs_core_simple(domain_cds, domain_protein):
    """analyze D/E vs core codon scores without CodonFM"""
    scores = compute_codon_nll_simple(domain_cds)

    de_scores = []
    core_scores = []

    for pos, codon, score in scores:
        if pos < len(domain_protein):
            aa = domain_protein[pos]
            if aa in 'DE':
                de_scores.append(score)
            else:
                core_scores.append(score)

    return de_scores, core_scores


def main():
    """main analysis pipeline"""
    CODONFM_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("EXPANDED CODONFM ANALYSIS")
    print("=" * 70)

    # load domain data
    print("\nLoading domain data...")
    domains = []
    with open(DATA_DIR / "domain_sequences.jsonl") as f:
        for line in f:
            d = json.loads(line)
            domains.append(d)

    # filter to human proteins with RRM_1 (largest RBD family)
    # human accessions typically start with P, Q, O, or are 6 chars
    human_rbd = [d for d in domains
                 if d['pfam_name'] in ['RRM_1', 'KH_1', 'DEAD', 'LSM', 'PUF']
                 and (d['accession'].startswith('P') or
                      d['accession'].startswith('Q') or
                      d['accession'].startswith('O'))
                 and len(d['accession']) == 6]

    # deduplicate by accession (keep first domain per protein)
    seen = set()
    unique_domains = []
    for d in human_rbd:
        key = (d['accession'], d['start'], d['end'])
        if key not in seen:
            seen.add(key)
            unique_domains.append(d)

    print(f"Human RBD domains: {len(unique_domains)}")

    # sample for analysis
    sample = unique_domains[:200]  # try 200
    print(f"Analyzing {len(sample)} domains...")

    # fetch CDS data
    valid_domains = []
    failed = 0

    print("\nFetching and validating CDS sequences...")

    for i, dom in enumerate(sample):
        if (i + 1) % 20 == 0:
            print(f"  Progress: {i+1}/{len(sample)}")

        # fetch UniProt data
        uniprot_data = fetch_cds_from_uniprot(dom['accession'])
        if not uniprot_data:
            failed += 1
            continue

        # try to get CDS from EMBL refs
        cds_found = False
        for ref in uniprot_data['cds_refs'][:3]:  # try first 3
            embl_id = ref.get('nucleotide_id', '').split('.')[0]
            if embl_id:
                full_cds = fetch_cds_sequence(embl_id)
                if full_cds:
                    # validate
                    valid, msg = validate_cds(full_cds, uniprot_data['protein_sequence'])
                    if valid:
                        # extract domain CDS
                        domain_cds = extract_domain_cds(full_cds, dom['start'], dom['end'])
                        if domain_cds and len(domain_cds) == len(dom['sequence']) * 3:
                            # verify translation matches
                            translated = translate(domain_cds)
                            if '*' not in translated and translated == dom['sequence']:
                                valid_domains.append({
                                    **dom,
                                    'domain_cds': domain_cds,
                                    'full_cds': full_cds,
                                    'cds_source': ref['database']
                                })
                                cds_found = True
                                break

        if not cds_found:
            failed += 1

        time.sleep(0.1)  # rate limiting

    print(f"\nValid domains with CDS: {len(valid_domains)}")
    print(f"Failed: {failed}")

    if len(valid_domains) < 20:
        print("ERROR: Not enough valid domains for analysis")
        return

    # analyze D/E vs core codon optimization
    print("\nAnalyzing D/E vs core codon optimization...")

    all_de_scores = []
    all_core_scores = []
    per_domain_results = []

    for dom in valid_domains:
        de_scores, core_scores = analyze_de_vs_core_simple(
            dom['domain_cds'], dom['sequence']
        )

        if de_scores and core_scores:
            all_de_scores.extend(de_scores)
            all_core_scores.extend(core_scores)

            per_domain_results.append({
                'accession': dom['accession'],
                'pfam': dom['pfam_name'],
                'de_mean': sum(de_scores) / len(de_scores),
                'core_mean': sum(core_scores) / len(core_scores),
                'n_de': len(de_scores),
                'n_core': len(core_scores),
            })

    # compute statistics
    import numpy as np
    from scipy import stats

    de_mean = np.mean(all_de_scores)
    core_mean = np.mean(all_core_scores)
    diff = de_mean - core_mean

    # statistical test
    t_stat, p_val = stats.ttest_ind(all_de_scores, all_core_scores)

    # effect size (Cohen's d)
    pooled_std = np.sqrt((np.std(all_de_scores)**2 + np.std(all_core_scores)**2) / 2)
    cohens_d = diff / pooled_std if pooled_std > 0 else 0

    print("\n" + "=" * 70)
    print("RESULTS: D/E vs CORE CODON OPTIMIZATION (CAI-based)")
    print("=" * 70)

    print(f"\nDomains analyzed: {len(per_domain_results)}")
    print(f"D/E codons: {len(all_de_scores)}")
    print(f"Core codons: {len(all_core_scores)}")

    print(f"\nMean codon score (negative frequency):")
    print(f"  D/E positions:  {de_mean:.4f} ± {np.std(all_de_scores):.4f}")
    print(f"  Core positions: {core_mean:.4f} ± {np.std(all_core_scores):.4f}")
    print(f"  Difference:     {diff:+.4f}")

    print(f"\nStatistical test:")
    print(f"  t-statistic: {t_stat:.4f}")
    print(f"  p-value: {p_val:.2e}")
    print(f"  Cohen's d: {cohens_d:.4f}")

    # interpretation
    if diff < 0:
        interpretation = "D/E codons use MORE COMMON codons (higher optimization)"
    else:
        interpretation = "D/E codons use RARER codons (lower optimization)"

    print(f"\nInterpretation: {interpretation}")

    # per-domain summary
    n_de_higher = sum(1 for r in per_domain_results if r['de_mean'] > r['core_mean'])
    print(f"\nDomains where D/E codons less optimal: {n_de_higher}/{len(per_domain_results)} ({100*n_de_higher/len(per_domain_results):.1f}%)")

    # save results
    output = {
        'n_domains': len(per_domain_results),
        'n_de_codons': len(all_de_scores),
        'n_core_codons': len(all_core_scores),
        'de_mean': float(de_mean),
        'core_mean': float(core_mean),
        'difference': float(diff),
        'p_value': float(p_val),
        'cohens_d': float(cohens_d),
        'per_domain': per_domain_results
    }

    with open(CODONFM_DIR / "expanded_analysis.json", "w") as f:
        json.dump(output, f, indent=2)

    # write summary
    with open(CODONFM_DIR / "expanded_analysis_summary.txt", "w") as f:
        f.write("=" * 70 + "\n")
        f.write("EXPANDED CODON ANALYSIS: D/E vs CORE POSITIONS\n")
        f.write("=" * 70 + "\n\n")

        f.write(f"Method: Codon Adaptation Index (CAI)-based scoring\n")
        f.write(f"Uses human codon usage frequencies from Kazusa database\n\n")

        f.write(f"Sample size:\n")
        f.write(f"  Domains: {len(per_domain_results)}\n")
        f.write(f"  D/E codons: {len(all_de_scores)}\n")
        f.write(f"  Core codons: {len(all_core_scores)}\n\n")

        f.write(f"Results:\n")
        f.write(f"  D/E mean score:  {de_mean:.4f} (higher = rarer codons)\n")
        f.write(f"  Core mean score: {core_mean:.4f}\n")
        f.write(f"  Difference: {diff:+.4f}\n\n")

        f.write(f"Statistics:\n")
        f.write(f"  p-value: {p_val:.2e}\n")
        f.write(f"  Cohen's d: {cohens_d:.4f}\n\n")

        f.write(f"Interpretation: {interpretation}\n\n")

        if abs(cohens_d) < 0.2:
            effect = "negligible"
        elif abs(cohens_d) < 0.5:
            effect = "small"
        elif abs(cohens_d) < 0.8:
            effect = "medium"
        else:
            effect = "large"

        f.write(f"Effect size: {effect}\n")

    print(f"\nResults saved to {CODONFM_DIR}/")


if __name__ == "__main__":
    main()
