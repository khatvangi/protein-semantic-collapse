#!/usr/bin/env python3
"""
PTM Enrichment Analysis in D/E Regions of RNA-Binding Proteins

Analyzes whether D/E (Asp/Glu) positions and D/E-rich segments in RBD domains
are enriched for post-translational modifications (phosphorylation, acetylation,
ubiquitination) compared to non-D/E positions.

Data source: UniProt REST API (freely accessible)
"""

import json
import requests
import time
import sys
from collections import defaultdict
from pathlib import Path
import re
from scipy import stats
import numpy as np

# configuration
DOMAIN_FILE = "/storage/kiran-stuff/protein-semantic-collapse/data/domain_sequences.jsonl"
OUTPUT_FILE = "/storage/kiran-stuff/protein-semantic-collapse/results/ptm_de_analysis.txt"
CACHE_DIR = Path("/storage/kiran-stuff/protein-semantic-collapse/data/ptm_cache")

# ptm categories of interest
PTM_CATEGORIES = {
    'phosphorylation': ['Phosphoserine', 'Phosphothreonine', 'Phosphotyrosine', 'Phospho'],
    'acetylation': ['N6-acetyllysine', 'acetyl', 'Acetyl'],
    'ubiquitination': ['ubiquitin', 'Ubiquitin', 'GlyGly'],
    'methylation': ['methyl', 'Methyl', 'dimethyl', 'trimethyl'],
}

def load_domain_data():
    """load domain sequences and group by accession"""
    domains_by_acc = defaultdict(list)
    with open(DOMAIN_FILE) as f:
        for line in f:
            data = json.loads(line)
            domains_by_acc[data['accession']].append(data)
    return domains_by_acc

def fetch_uniprot_ptm(accession, cache_dir=CACHE_DIR):
    """fetch PTM annotations from UniProt for a single accession"""
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_file = cache_dir / f"{accession}.json"

    # use cache if available
    if cache_file.exists():
        with open(cache_file) as f:
            return json.load(f)

    url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
    try:
        resp = requests.get(url, timeout=30)
        if resp.status_code == 200:
            data = resp.json()
            # extract PTM features and sequence
            ptm_data = {
                'accession': accession,
                'sequence': data.get('sequence', {}).get('value', ''),
                'ptms': []
            }
            for feature in data.get('features', []):
                ftype = feature.get('type', '')
                if ftype in ['Modified residue', 'Cross-link', 'Lipidation', 'Glycosylation']:
                    loc = feature.get('location', {})
                    start = loc.get('start', {}).get('value')
                    end = loc.get('end', {}).get('value')
                    desc = feature.get('description', '')
                    if start and end:
                        ptm_data['ptms'].append({
                            'type': ftype,
                            'position': start,  # for single-residue PTMs start==end
                            'description': desc
                        })

            # cache result
            with open(cache_file, 'w') as f:
                json.dump(ptm_data, f)
            return ptm_data
        else:
            return None
    except Exception as e:
        print(f"  Error fetching {accession}: {e}", file=sys.stderr)
        return None

def categorize_ptm(description):
    """categorize PTM by type based on description"""
    desc_lower = description.lower()
    categories = []
    for cat, keywords in PTM_CATEGORIES.items():
        for kw in keywords:
            if kw.lower() in desc_lower:
                categories.append(cat)
                break
    return categories if categories else ['other']

def find_de_rich_segments(sequence, min_length=3):
    """find D/E-rich segments (>=min_length consecutive D/E)"""
    segments = []
    i = 0
    while i < len(sequence):
        if sequence[i] in 'DE':
            start = i
            while i < len(sequence) and sequence[i] in 'DE':
                i += 1
            if i - start >= min_length:
                segments.append((start, i))  # 0-indexed, exclusive end
        else:
            i += 1
    return segments

def analyze_protein(accession, domains, ptm_data):
    """analyze PTM enrichment for a single protein"""
    if ptm_data is None or not ptm_data.get('sequence'):
        return None

    full_seq = ptm_data['sequence']
    results = {
        'accession': accession,
        'domain_ptms': defaultdict(lambda: {'de': [], 'non_de': [], 'de_rich': [], 'domain_core': []}),
        'domain_info': []
    }

    # process each domain
    for domain in domains:
        dom_start = domain['start'] - 1  # convert to 0-indexed
        dom_end = domain['end']  # keep as is for slicing
        dom_seq = domain['sequence']
        dom_name = domain['pfam_name']

        # validate domain position matches sequence
        if dom_end > len(full_seq):
            continue

        # identify D/E positions within domain (relative to full protein)
        de_positions = set()
        non_de_positions = set()
        for i, aa in enumerate(dom_seq):
            global_pos = dom_start + i + 1  # 1-indexed position in full protein
            if aa in 'DE':
                de_positions.add(global_pos)
            else:
                non_de_positions.add(global_pos)

        # identify D/E-rich segments
        de_rich_segments = find_de_rich_segments(dom_seq)
        de_rich_positions = set()
        for seg_start, seg_end in de_rich_segments:
            for j in range(seg_start, seg_end):
                de_rich_positions.add(dom_start + j + 1)

        # domain core = non-DE positions outside DE-rich segments
        domain_core_positions = non_de_positions - de_rich_positions

        results['domain_info'].append({
            'pfam': domain['pfam'],
            'pfam_name': dom_name,
            'start': domain['start'],
            'end': domain['end'],
            'length': domain['domain_length'],
            'de_count': len(de_positions),
            'de_rich_segments': len(de_rich_segments),
            'de_fraction': len(de_positions) / domain['domain_length']
        })

        # map PTMs to regions
        for ptm in ptm_data['ptms']:
            pos = ptm['position']
            if dom_start + 1 <= pos <= dom_end:
                cats = categorize_ptm(ptm['description'])
                for cat in cats:
                    if pos in de_positions:
                        results['domain_ptms'][cat]['de'].append(pos)
                    else:
                        results['domain_ptms'][cat]['non_de'].append(pos)

                    if pos in de_rich_positions:
                        results['domain_ptms'][cat]['de_rich'].append(pos)
                    elif pos in domain_core_positions:
                        results['domain_ptms'][cat]['domain_core'].append(pos)

        # store position counts for later statistical analysis
        results[f'domain_{domain["start"]}_{domain["end"]}'] = {
            'de_positions': len(de_positions),
            'non_de_positions': len(non_de_positions),
            'de_rich_positions': len(de_rich_positions),
            'domain_core_positions': len(domain_core_positions)
        }

    return results

def aggregate_results(all_results):
    """aggregate results across all proteins"""
    agg = {
        'total_proteins': 0,
        'total_domains': 0,
        'ptm_counts': defaultdict(lambda: {'de': 0, 'non_de': 0, 'de_rich': 0, 'domain_core': 0}),
        'position_counts': {'de': 0, 'non_de': 0, 'de_rich': 0, 'domain_core': 0},
        'proteins_with_ptms': 0,
        'de_fraction_distribution': []
    }

    for result in all_results:
        if result is None:
            continue
        agg['total_proteins'] += 1
        agg['total_domains'] += len(result['domain_info'])

        has_ptm = False
        for cat, regions in result['domain_ptms'].items():
            for region, ptms in regions.items():
                agg['ptm_counts'][cat][region] += len(ptms)
                if len(ptms) > 0:
                    has_ptm = True

        if has_ptm:
            agg['proteins_with_ptms'] += 1

        # aggregate position counts
        for key in result:
            if key.startswith('domain_') and key != 'domain_ptms' and key != 'domain_info':
                counts = result[key]
                agg['position_counts']['de'] += counts['de_positions']
                agg['position_counts']['non_de'] += counts['non_de_positions']
                agg['position_counts']['de_rich'] += counts['de_rich_positions']
                agg['position_counts']['domain_core'] += counts['domain_core_positions']

        for dom_info in result['domain_info']:
            agg['de_fraction_distribution'].append(dom_info['de_fraction'])

    return agg

def calculate_enrichment(ptm_count, total_positions, baseline_count, baseline_positions):
    """calculate enrichment ratio and p-value using Fisher's exact test"""
    if total_positions == 0 or baseline_positions == 0:
        return None, None, None

    # contingency table:
    # [[PTM in region, non-PTM in region], [PTM in baseline, non-PTM in baseline]]
    table = [
        [ptm_count, total_positions - ptm_count],
        [baseline_count, baseline_positions - baseline_count]
    ]

    # observed rate vs expected rate
    rate_region = ptm_count / total_positions if total_positions > 0 else 0
    rate_baseline = baseline_count / baseline_positions if baseline_positions > 0 else 0

    enrichment = (rate_region / rate_baseline) if rate_baseline > 0 else float('inf')

    # fisher's exact test
    try:
        odds_ratio, pvalue = stats.fisher_exact(table)
    except:
        odds_ratio, pvalue = None, None

    return enrichment, odds_ratio, pvalue

def write_report(agg, output_file):
    """write comprehensive analysis report"""
    with open(output_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("PTM ENRICHMENT ANALYSIS IN D/E REGIONS OF RNA-BINDING DOMAINS\n")
        f.write("=" * 80 + "\n\n")

        # overview
        f.write("OVERVIEW\n")
        f.write("-" * 40 + "\n")
        f.write(f"Total proteins analyzed: {agg['total_proteins']}\n")
        f.write(f"Total domains analyzed: {agg['total_domains']}\n")
        f.write(f"Proteins with at least one PTM: {agg['proteins_with_ptms']}\n\n")

        # position statistics
        f.write("POSITION STATISTICS (within domains)\n")
        f.write("-" * 40 + "\n")
        total_pos = agg['position_counts']['de'] + agg['position_counts']['non_de']
        f.write(f"Total positions in domains: {total_pos}\n")
        f.write(f"D/E positions: {agg['position_counts']['de']} ({100*agg['position_counts']['de']/total_pos:.1f}%)\n")
        f.write(f"Non-D/E positions: {agg['position_counts']['non_de']} ({100*agg['position_counts']['non_de']/total_pos:.1f}%)\n")
        f.write(f"D/E-rich segment positions (>=3 consecutive D/E): {agg['position_counts']['de_rich']}\n")
        f.write(f"Domain core positions (non-D/E outside D/E-rich): {agg['position_counts']['domain_core']}\n\n")

        # D/E fraction distribution
        de_fracs = agg['de_fraction_distribution']
        if de_fracs:
            f.write("D/E FRACTION PER DOMAIN\n")
            f.write("-" * 40 + "\n")
            f.write(f"Mean D/E fraction: {np.mean(de_fracs):.3f}\n")
            f.write(f"Median D/E fraction: {np.median(de_fracs):.3f}\n")
            f.write(f"Std D/E fraction: {np.std(de_fracs):.3f}\n")
            f.write(f"Min/Max: {np.min(de_fracs):.3f} / {np.max(de_fracs):.3f}\n\n")

        # PTM counts by category
        f.write("PTM COUNTS BY CATEGORY AND REGION\n")
        f.write("-" * 40 + "\n")
        f.write(f"{'Category':<20} {'D/E':>8} {'Non-D/E':>10} {'D/E-rich':>10} {'Core':>8}\n")
        f.write("-" * 60 + "\n")

        all_de = 0
        all_non_de = 0
        all_de_rich = 0
        all_core = 0

        for cat in ['phosphorylation', 'acetylation', 'ubiquitination', 'methylation', 'other']:
            if cat in agg['ptm_counts']:
                counts = agg['ptm_counts'][cat]
                f.write(f"{cat:<20} {counts['de']:>8} {counts['non_de']:>10} {counts['de_rich']:>10} {counts['domain_core']:>8}\n")
                all_de += counts['de']
                all_non_de += counts['non_de']
                all_de_rich += counts['de_rich']
                all_core += counts['domain_core']

        f.write("-" * 60 + "\n")
        f.write(f"{'TOTAL':<20} {all_de:>8} {all_non_de:>10} {all_de_rich:>10} {all_core:>8}\n\n")

        # enrichment analysis
        f.write("ENRICHMENT ANALYSIS\n")
        f.write("=" * 80 + "\n\n")

        # 1. D/E vs non-D/E positions
        f.write("1. PTM DENSITY: D/E POSITIONS vs NON-D/E POSITIONS\n")
        f.write("-" * 60 + "\n")
        de_pos = agg['position_counts']['de']
        non_de_pos = agg['position_counts']['non_de']

        f.write(f"{'Category':<20} {'D/E rate':>12} {'Non-D/E rate':>14} {'Enrichment':>12} {'p-value':>12}\n")
        f.write("-" * 70 + "\n")

        for cat in ['phosphorylation', 'acetylation', 'ubiquitination', 'methylation', 'other']:
            if cat in agg['ptm_counts']:
                counts = agg['ptm_counts'][cat]
                enrichment, odds_ratio, pval = calculate_enrichment(
                    counts['de'], de_pos, counts['non_de'], non_de_pos
                )

                de_rate = counts['de'] / de_pos * 1000 if de_pos > 0 else 0  # per 1000 positions
                non_de_rate = counts['non_de'] / non_de_pos * 1000 if non_de_pos > 0 else 0

                enr_str = f"{enrichment:.2f}x" if enrichment and enrichment != float('inf') else "N/A"
                pval_str = f"{pval:.2e}" if pval else "N/A"

                f.write(f"{cat:<20} {de_rate:>10.2f}/k {non_de_rate:>12.2f}/k {enr_str:>12} {pval_str:>12}\n")

        # total across all PTM types
        enrichment_total, odds_total, pval_total = calculate_enrichment(
            all_de, de_pos, all_non_de, non_de_pos
        )
        de_rate_total = all_de / de_pos * 1000 if de_pos > 0 else 0
        non_de_rate_total = all_non_de / non_de_pos * 1000 if non_de_pos > 0 else 0
        enr_str = f"{enrichment_total:.2f}x" if enrichment_total and enrichment_total != float('inf') else "N/A"
        pval_str = f"{pval_total:.2e}" if pval_total else "N/A"

        f.write("-" * 70 + "\n")
        f.write(f"{'ALL PTMs':<20} {de_rate_total:>10.2f}/k {non_de_rate_total:>12.2f}/k {enr_str:>12} {pval_str:>12}\n\n")

        # 2. D/E-rich segments vs domain core
        f.write("2. PTM DENSITY: D/E-RICH SEGMENTS (>=3 D/E) vs DOMAIN CORE\n")
        f.write("-" * 60 + "\n")
        de_rich_pos = agg['position_counts']['de_rich']
        core_pos = agg['position_counts']['domain_core']

        f.write(f"{'Category':<20} {'D/E-rich rate':>14} {'Core rate':>12} {'Enrichment':>12} {'p-value':>12}\n")
        f.write("-" * 70 + "\n")

        for cat in ['phosphorylation', 'acetylation', 'ubiquitination', 'methylation', 'other']:
            if cat in agg['ptm_counts']:
                counts = agg['ptm_counts'][cat]
                enrichment, odds_ratio, pval = calculate_enrichment(
                    counts['de_rich'], de_rich_pos, counts['domain_core'], core_pos
                )

                de_rich_rate = counts['de_rich'] / de_rich_pos * 1000 if de_rich_pos > 0 else 0
                core_rate = counts['domain_core'] / core_pos * 1000 if core_pos > 0 else 0

                enr_str = f"{enrichment:.2f}x" if enrichment and enrichment != float('inf') else "N/A"
                pval_str = f"{pval:.2e}" if pval else "N/A"

                f.write(f"{cat:<20} {de_rich_rate:>12.2f}/k {core_rate:>10.2f}/k {enr_str:>12} {pval_str:>12}\n")

        f.write("\n")

        # interpretation
        f.write("INTERPRETATION\n")
        f.write("=" * 80 + "\n")
        f.write("""
The analysis compares PTM density between:
1. D/E positions (Asp/Glu) vs non-D/E positions within RBD domains
2. D/E-rich segments (>=3 consecutive D/E) vs domain core regions

Key metrics:
- Rate: PTMs per 1000 positions (‰)
- Enrichment: Ratio of D/E rate to non-D/E rate (>1 = enriched in D/E)
- p-value: Fisher's exact test (significant if <0.05)

Note: D and E residues are NOT typical phosphorylation targets (S, T, Y are).
However, D/E-rich regions may be regulatory hotspots for OTHER PTMs or serve
as PTM-adjacent regulatory sequences.

If D/E-rich regions show PTM enrichment, this supports the hypothesis that
D/E segments are regulatory hotspots in RNA-binding proteins.
""")

        # data source
        f.write("\nDATA SOURCE\n")
        f.write("-" * 40 + "\n")
        f.write("PTM annotations: UniProt REST API (https://rest.uniprot.org)\n")
        f.write("Domain data: /storage/kiran-stuff/protein-semantic-collapse/data/domain_sequences.jsonl\n")
        f.write(f"Analysis date: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")

def main():
    print("Loading domain data...")
    domains_by_acc = load_domain_data()
    print(f"Found {len(domains_by_acc)} unique proteins with {sum(len(v) for v in domains_by_acc.values())} domains")

    # filter to human/well-annotated proteins for better PTM coverage
    # UniProt accessions starting with certain patterns indicate species
    accessions = list(domains_by_acc.keys())

    print(f"\nFetching PTM data from UniProt for {len(accessions)} proteins...")
    print("(Using cache when available)")

    all_results = []
    fetched = 0
    cached = 0
    failed = 0

    for i, acc in enumerate(accessions):
        if (i + 1) % 100 == 0:
            print(f"  Progress: {i+1}/{len(accessions)} (fetched: {fetched}, cached: {cached}, failed: {failed})")

        cache_file = CACHE_DIR / f"{acc}.json"
        was_cached = cache_file.exists()

        ptm_data = fetch_uniprot_ptm(acc)

        if ptm_data:
            if was_cached:
                cached += 1
            else:
                fetched += 1
                time.sleep(0.1)  # rate limiting for API

            result = analyze_protein(acc, domains_by_acc[acc], ptm_data)
            all_results.append(result)
        else:
            failed += 1

    print(f"\nFetch complete: {fetched} new, {cached} cached, {failed} failed")

    print("\nAggregating results...")
    agg = aggregate_results(all_results)

    print("\nWriting report...")
    write_report(agg, OUTPUT_FILE)
    print(f"Report written to: {OUTPUT_FILE}")

    # print summary to stdout
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Proteins analyzed: {agg['total_proteins']}")
    print(f"Domains analyzed: {agg['total_domains']}")
    print(f"Proteins with PTMs: {agg['proteins_with_ptms']}")
    print(f"Total PTMs in D/E positions: {sum(c['de'] for c in agg['ptm_counts'].values())}")
    print(f"Total PTMs in non-D/E positions: {sum(c['non_de'] for c in agg['ptm_counts'].values())}")

if __name__ == "__main__":
    main()
