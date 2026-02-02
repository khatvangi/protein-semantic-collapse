#!/usr/bin/env python3
"""
Evolutionary Analysis of D/E Autoinhibition in RBD Families

Questions:
1. Which RBD families are ancient (all domains of life) vs derived (eukaryotes only)?
2. Is there a correlation between domain age and D/E content?
3. When did D/E autoinhibition emerge?

Approach:
1. Use Pfam/InterPro taxonomic distribution data
2. Classify families by phylogenetic breadth
3. Map D/E content onto evolutionary timeline
4. Perform phylogenetic analysis of representative sequences

Data sources:
- Our domain data (sequences, D/E content)
- Pfam taxonomic distributions (API)
- NCBI taxonomy for species classification
"""

import json
import requests
from collections import defaultdict, Counter
from pathlib import Path
import numpy as np
from scipy import stats
import time

RESULTS_DIR = Path("/storage/kiran-stuff/protein-semantic-collapse/results")

# known evolutionary origins of RBD families (from literature/Pfam)
# categories:
#   'LUCA' = Last Universal Common Ancestor (all 3 domains)
#   'bacterial' = bacteria + archaea (no eukaryotes or rare)
#   'eukaryotic' = eukaryotes only
#   'metazoan' = animals only
FAMILY_ORIGINS = {
    # ancient (LUCA or bacterial origin)
    'KH_1': 'LUCA',           # found in bacteria, archaea, eukaryotes
    'S1': 'LUCA',             # ribosomal, all domains
    'CSD': 'bacterial',       # cold-shock, bacterial origin but spread
    'OB_NTP_bind': 'LUCA',    # ancient nucleotide binding
    'GTP_EFTU': 'LUCA',       # translation factor, universal
    'Ribosomal_S1': 'LUCA',   # ribosomal

    # eukaryotic innovations
    'RRM_1': 'eukaryotic',    # major eukaryotic expansion
    'LSM': 'LUCA',            # Sm/Lsm - actually ancient but expanded in eukaryotes
    'DEAD': 'eukaryotic',     # DEAD-box helicases, eukaryotic expansion
    'Helicase_C': 'eukaryotic',
    'dsRBD': 'eukaryotic',    # double-stranded RNA binding
    'PUF': 'eukaryotic',      # Pumilio, eukaryotic
    'PAZ': 'eukaryotic',      # RNA silencing, eukaryotic

    # zinc fingers - various origins
    'zf-CCHC': 'eukaryotic',  # retroviral/eukaryotic
    'zf-CCCH': 'eukaryotic',  # CCCH zinc fingers
    'zf-C2H2': 'eukaryotic',  # classic zinc finger, eukaryotic

    # other
    'WD40': 'LUCA',           # ancient scaffold
    'Pkinase': 'LUCA',        # protein kinases, ancient
}

def load_domains():
    """load domain data"""
    with open(RESULTS_DIR / "rbp_domains.json") as f:
        return json.load(f)

def calculate_de_content(seq):
    """calculate D/E percentage"""
    if not seq:
        return 0
    return sum(1 for aa in seq if aa in 'DE') / len(seq) * 100

def calculate_kr_content(seq):
    """calculate K/R (basic) percentage"""
    if not seq:
        return 0
    return sum(1 for aa in seq if aa in 'KR') / len(seq) * 100

def get_pfam_species_count(pfam_acc):
    """
    query Pfam API for species distribution
    returns count of species with this domain
    """
    # use InterPro API (Pfam is now part of InterPro)
    url = f"https://www.ebi.ac.uk/interpro/api/entry/pfam/{pfam_acc}/"

    try:
        response = requests.get(url, timeout=30)
        if response.status_code == 200:
            data = response.json()
            # get taxonomic distribution if available
            counters = data.get('metadata', {}).get('counters', {})
            return {
                'proteins': counters.get('proteins', 0),
                'species': counters.get('taxa', 0),
                'structures': counters.get('structures', 0)
            }
    except Exception as e:
        print(f"  Warning: Could not fetch data for {pfam_acc}: {e}")

    return None

def analyze_family_de(domains, family_name):
    """analyze D/E content for a family"""
    family_domains = [d for d in domains if d['pfam_name'] == family_name]

    if not family_domains:
        return None

    de_contents = [calculate_de_content(d['sequence']) for d in family_domains]
    kr_contents = [calculate_kr_content(d['sequence']) for d in family_domains]

    return {
        'family': family_name,
        'n_domains': len(family_domains),
        'mean_de': np.mean(de_contents),
        'std_de': np.std(de_contents),
        'mean_kr': np.mean(kr_contents),
        'de_kr_ratio': np.mean(de_contents) / np.mean(kr_contents) if np.mean(kr_contents) > 0 else 0
    }

def get_representative_sequences(domains, family_name, n_seqs=10):
    """get representative sequences for phylogenetic analysis"""
    family_domains = [d for d in domains if d['pfam_name'] == family_name]

    if len(family_domains) <= n_seqs:
        return family_domains

    # sample to get diversity - take from different D/E content ranges
    de_contents = [(d, calculate_de_content(d['sequence'])) for d in family_domains]
    de_contents.sort(key=lambda x: x[1])

    # take evenly spaced samples
    indices = np.linspace(0, len(de_contents) - 1, n_seqs, dtype=int)
    return [de_contents[i][0] for i in indices]

def write_fasta(sequences, output_file, family_name):
    """write sequences to FASTA file"""
    with open(output_file, 'w') as f:
        for i, seq in enumerate(sequences):
            acc = seq.get('accession', f'seq{i}')
            de_content = calculate_de_content(seq['sequence'])
            f.write(f">{acc}_{family_name}_DE{de_content:.1f}\n")
            f.write(f"{seq['sequence']}\n")

def main():
    print("=" * 70)
    print("Evolutionary Analysis of D/E Autoinhibition")
    print("=" * 70)

    # load data
    domains = load_domains()
    print(f"\nLoaded {len(domains)} domains")

    # get unique families and their sizes
    family_counts = Counter(d['pfam_name'] for d in domains)
    print(f"Total families: {len(family_counts)}")

    # analyze D/E content by family
    print("\n" + "-" * 70)
    print("D/E Content by Evolutionary Origin")
    print("-" * 70)

    family_analyses = {}
    for family in family_counts:
        analysis = analyze_family_de(domains, family)
        if analysis:
            family_analyses[family] = analysis

    # group by evolutionary origin
    origin_groups = defaultdict(list)
    unclassified = []

    for family, analysis in family_analyses.items():
        if family in FAMILY_ORIGINS:
            origin = FAMILY_ORIGINS[family]
            origin_groups[origin].append(analysis)
        else:
            unclassified.append(analysis)

    # summarize by origin
    print(f"\n{'Origin':<15} {'n_families':>10} {'n_domains':>10} {'Mean D/E':>10} {'Std D/E':>10}")
    print("-" * 55)

    origin_de_means = {}
    for origin in ['LUCA', 'bacterial', 'eukaryotic']:
        if origin in origin_groups:
            analyses = origin_groups[origin]
            n_fam = len(analyses)
            n_dom = sum(a['n_domains'] for a in analyses)
            # weighted mean by domain count
            total_de = sum(a['mean_de'] * a['n_domains'] for a in analyses)
            mean_de = total_de / n_dom if n_dom > 0 else 0
            # calculate std across families
            de_values = [a['mean_de'] for a in analyses]
            std_de = np.std(de_values) if len(de_values) > 1 else 0

            origin_de_means[origin] = de_values
            print(f"{origin:<15} {n_fam:>10} {n_dom:>10} {mean_de:>10.2f}% {std_de:>10.2f}")

    # statistical test: LUCA vs eukaryotic
    print("\n" + "-" * 70)
    print("Statistical Comparison: Ancient (LUCA) vs Eukaryotic")
    print("-" * 70)

    if 'LUCA' in origin_de_means and 'eukaryotic' in origin_de_means:
        luca_de = origin_de_means['LUCA']
        euk_de = origin_de_means['eukaryotic']

        stat, pval = stats.mannwhitneyu(luca_de, euk_de, alternative='two-sided')
        print(f"\nLUCA families mean D/E: {np.mean(luca_de):.2f}% (n={len(luca_de)})")
        print(f"Eukaryotic families mean D/E: {np.mean(euk_de):.2f}% (n={len(euk_de)})")
        print(f"Difference: {np.mean(euk_de) - np.mean(luca_de):+.2f}%")
        print(f"Mann-Whitney U p-value: {pval:.4f}")

        if pval < 0.05:
            print("\n*** SIGNIFICANT: Eukaryotic families have different D/E content than ancient families ***")

    # detailed family listing
    print("\n" + "-" * 70)
    print("Individual Family Analysis (sorted by D/E content)")
    print("-" * 70)

    sorted_families = sorted(family_analyses.values(), key=lambda x: x['mean_de'], reverse=True)

    print(f"\n{'Family':<20} {'Origin':<12} {'n_domains':>10} {'D/E %':>8} {'K/R %':>8} {'D/E:K/R':>8}")
    print("-" * 70)

    for analysis in sorted_families[:25]:  # top 25
        family = analysis['family']
        origin = FAMILY_ORIGINS.get(family, 'unknown')
        print(f"{family:<20} {origin:<12} {analysis['n_domains']:>10} {analysis['mean_de']:>8.2f} {analysis['mean_kr']:>8.2f} {analysis['de_kr_ratio']:>8.2f}")

    # key finding: which ancient families DO have high D/E?
    print("\n" + "=" * 70)
    print("KEY EVOLUTIONARY FINDINGS")
    print("=" * 70)

    # check LSM specifically - ancient but high D/E
    lsm_analysis = family_analyses.get('LSM')
    if lsm_analysis:
        print(f"\n1. LSM family (ancient, LUCA origin):")
        print(f"   D/E content: {lsm_analysis['mean_de']:.2f}%")
        print(f"   This is ANCIENT yet uses D/E autoinhibition!")

    # check KH specifically - ancient but low D/E
    kh_analysis = family_analyses.get('KH_1')
    if kh_analysis:
        print(f"\n2. KH_1 family (ancient, LUCA origin):")
        print(f"   D/E content: {kh_analysis['mean_de']:.2f}%")
        print(f"   Ancient but does NOT use D/E autoinhibition")

    # RRM - eukaryotic with high D/E
    rrm_analysis = family_analyses.get('RRM_1')
    if rrm_analysis:
        print(f"\n3. RRM_1 family (eukaryotic innovation):")
        print(f"   D/E content: {rrm_analysis['mean_de']:.2f}%")
        print(f"   Eukaryotic expansion with D/E autoinhibition")

    # prepare sequences for phylogenetic analysis
    print("\n" + "-" * 70)
    print("Preparing sequences for phylogenetic analysis...")
    print("-" * 70)

    phylo_dir = RESULTS_DIR / "phylogenetics"
    phylo_dir.mkdir(exist_ok=True)

    # write representative sequences for key families
    key_families = ['RRM_1', 'KH_1', 'LSM', 'DEAD', 'dsRBD', 'S1']
    available_families = set(family_analyses.keys())

    for family in key_families:
        if family in available_families:
            seqs = get_representative_sequences(domains, family, n_seqs=20)
            fasta_file = phylo_dir / f"{family}_representatives.fasta"
            write_fasta(seqs, fasta_file, family)
            print(f"  Written {len(seqs)} sequences for {family} to {fasta_file.name}")

    # save full results
    output = {
        'family_analyses': {k: v for k, v in family_analyses.items()},
        'origin_groups': {k: [a['family'] for a in v] for k, v in origin_groups.items()},
        'origin_de_means': {k: {'mean': float(np.mean(v)), 'std': float(np.std(v)), 'n': len(v)}
                           for k, v in origin_de_means.items()},
        'unclassified_families': [a['family'] for a in unclassified]
    }

    with open(RESULTS_DIR / "evolutionary_analysis.json", 'w') as f:
        json.dump(output, f, indent=2, default=str)

    print(f"\nResults saved to {RESULTS_DIR / 'evolutionary_analysis.json'}")
    print(f"FASTA files for phylogenetic analysis in {phylo_dir}/")

    return family_analyses

if __name__ == "__main__":
    main()
