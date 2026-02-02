#!/usr/bin/env python3
"""
Deep Functional Analysis: D/E in Dynamic vs Stable RBPs

Four analyses:
1. mRNA decay proteins - do they have high D/E like splicing?
2. SR protein phosphorylation - D/E + phospho-S/R charge regulation
3. Stress granule vs ribosomal proteins - dynamic vs stable
4. KH regulatory mechanisms - what do non-D/E families use?
"""

import json
import requests
import numpy as np
from scipy import stats
from collections import Counter, defaultdict
import time
from pathlib import Path

RESULTS_DIR = Path("/storage/kiran-stuff/protein-semantic-collapse/results")

# load domain data
with open(RESULTS_DIR / "rbp_domains.json") as f:
    domains = json.load(f)

def calculate_de(seq):
    return sum(1 for aa in seq if aa in 'DE') / len(seq) * 100 if seq else 0

def calculate_kr(seq):
    return sum(1 for aa in seq if aa in 'KR') / len(seq) * 100 if seq else 0

def calculate_sr(seq):
    """serine + arginine content"""
    return sum(1 for aa in seq if aa in 'SR') / len(seq) * 100 if seq else 0

def get_protein_go_and_features(accession):
    """fetch GO terms and features from UniProt"""
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
    try:
        resp = requests.get(url, timeout=10)
        if resp.status_code == 200:
            data = resp.json()

            # GO terms
            go_terms = []
            for xref in data.get('uniProtKBCrossReferences', []):
                if xref.get('database') == 'GO':
                    props = {p['key']: p['value'] for p in xref.get('properties', [])}
                    go_terms.append(props.get('GoTerm', ''))

            # protein name
            name = data.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', '')

            # features (PTMs)
            features = data.get('features', [])
            phospho_sites = [f for f in features if f.get('type') == 'Modified residue'
                           and 'phospho' in f.get('description', '').lower()]

            # subcellular location
            comments = data.get('comments', [])
            locations = []
            for c in comments:
                if c.get('commentType') == 'SUBCELLULAR LOCATION':
                    for loc in c.get('subcellularLocations', []):
                        locations.append(loc.get('location', {}).get('value', ''))

            return {
                'go_terms': go_terms,
                'name': name,
                'phospho_count': len(phospho_sites),
                'locations': locations
            }
    except Exception as e:
        pass
    return {'go_terms': [], 'name': '', 'phospho_count': 0, 'locations': []}

# ============================================================
# ANALYSIS 1: mRNA Decay Proteins
# ============================================================
def analyze_mrna_decay():
    print("=" * 70)
    print("ANALYSIS 1: mRNA Decay Proteins")
    print("=" * 70)

    decay_keywords = ['decay', 'degradation', 'exosome', 'deadenylation',
                      'decapping', 'P-body', 'XRN', 'DCP', 'CCR4', 'NOT']
    splicing_keywords = ['splicing', 'spliceosome']

    # get sample of proteins
    all_accessions = list(set(d['accession'] for d in domains))[:100]

    decay_proteins = []
    splicing_proteins = []
    other_proteins = []

    print(f"\nClassifying {len(all_accessions)} proteins...")

    for i, acc in enumerate(all_accessions):
        info = get_protein_go_and_features(acc)
        combined = ' '.join(info['go_terms'] + [info['name']]).lower()

        # get D/E for this protein
        acc_doms = [d for d in domains if d['accession'] == acc]
        de_pcts = [calculate_de(d['sequence']) for d in acc_doms]
        mean_de = np.mean(de_pcts) if de_pcts else 0

        is_decay = any(kw.lower() in combined for kw in decay_keywords)
        is_splicing = any(kw.lower() in combined for kw in splicing_keywords)

        if is_decay:
            decay_proteins.append({'acc': acc, 'de': mean_de, 'name': info['name'][:50]})
        elif is_splicing:
            splicing_proteins.append({'acc': acc, 'de': mean_de, 'name': info['name'][:50]})
        else:
            other_proteins.append({'acc': acc, 'de': mean_de, 'name': info['name'][:50]})

        if (i+1) % 25 == 0:
            print(f"  {i+1}/{len(all_accessions)}")
        time.sleep(0.1)

    print(f"\nResults:")
    print(f"  Decay proteins: {len(decay_proteins)}")
    print(f"  Splicing proteins: {len(splicing_proteins)}")
    print(f"  Other proteins: {len(other_proteins)}")

    if decay_proteins:
        decay_de = [p['de'] for p in decay_proteins]
        print(f"\n  Decay D/E: {np.mean(decay_de):.2f}% ± {np.std(decay_de):.2f}")
    if splicing_proteins:
        splicing_de = [p['de'] for p in splicing_proteins]
        print(f"  Splicing D/E: {np.mean(splicing_de):.2f}% ± {np.std(splicing_de):.2f}")
    if other_proteins:
        other_de = [p['de'] for p in other_proteins]
        print(f"  Other D/E: {np.mean(other_de):.2f}% ± {np.std(other_de):.2f}")

    # statistical test
    if decay_proteins and other_proteins and len(decay_proteins) >= 3:
        stat, pval = stats.mannwhitneyu([p['de'] for p in decay_proteins],
                                        [p['de'] for p in other_proteins])
        print(f"\n  Decay vs Other: p = {pval:.4f}")
        if pval < 0.05:
            print("  *** DECAY PROTEINS HAVE DIFFERENT D/E ***")

    return {'decay': decay_proteins, 'splicing': splicing_proteins, 'other': other_proteins}

# ============================================================
# ANALYSIS 2: SR Protein Phosphorylation
# ============================================================
def analyze_sr_phosphorylation():
    print("\n" + "=" * 70)
    print("ANALYSIS 2: SR Protein Phosphorylation + D/E")
    print("=" * 70)

    # SR proteins are serine/arginine-rich splicing factors
    # they have high S+R content AND are heavily phosphorylated

    # get RRM proteins and check for SR-rich ones
    rrm_accessions = list(set(d['accession'] for d in domains if d['pfam_name'] == 'RRM_1'))[:50]

    sr_proteins = []
    non_sr_proteins = []

    print(f"\nAnalyzing {len(rrm_accessions)} RRM proteins for SR content and phosphorylation...")

    for i, acc in enumerate(rrm_accessions):
        info = get_protein_go_and_features(acc)

        # get sequence composition
        acc_doms = [d for d in domains if d['accession'] == acc]
        all_seq = ''.join(d['sequence'] for d in acc_doms)

        de_pct = calculate_de(all_seq)
        sr_pct = calculate_sr(all_seq)

        # SR proteins typically have >15% S+R content
        is_sr = sr_pct > 15 or 'SRSF' in info['name'].upper() or 'SR ' in info['name'].upper()

        protein_data = {
            'acc': acc,
            'name': info['name'][:40],
            'de': de_pct,
            'sr': sr_pct,
            'phospho': info['phospho_count']
        }

        if is_sr:
            sr_proteins.append(protein_data)
        else:
            non_sr_proteins.append(protein_data)

        if (i+1) % 25 == 0:
            print(f"  {i+1}/{len(rrm_accessions)}")
        time.sleep(0.1)

    print(f"\nResults:")
    print(f"  SR proteins: {len(sr_proteins)}")
    print(f"  Non-SR proteins: {len(non_sr_proteins)}")

    if sr_proteins:
        print(f"\n  SR proteins:")
        print(f"    D/E: {np.mean([p['de'] for p in sr_proteins]):.2f}%")
        print(f"    S+R: {np.mean([p['sr'] for p in sr_proteins]):.2f}%")
        print(f"    Phospho sites: {np.mean([p['phospho'] for p in sr_proteins]):.1f}")

    if non_sr_proteins:
        print(f"\n  Non-SR proteins:")
        print(f"    D/E: {np.mean([p['de'] for p in non_sr_proteins]):.2f}%")
        print(f"    S+R: {np.mean([p['sr'] for p in non_sr_proteins]):.2f}%")
        print(f"    Phospho sites: {np.mean([p['phospho'] for p in non_sr_proteins]):.1f}")

    # correlation between D/E and phosphorylation
    all_proteins = sr_proteins + non_sr_proteins
    if len(all_proteins) >= 10:
        de_vals = [p['de'] for p in all_proteins]
        phospho_vals = [p['phospho'] for p in all_proteins]

        corr, pval = stats.spearmanr(de_vals, phospho_vals)
        print(f"\n  Correlation (D/E vs Phospho): r = {corr:.3f}, p = {pval:.4f}")

        if pval < 0.05:
            print("  *** D/E CORRELATES WITH PHOSPHORYLATION ***")

    return {'sr': sr_proteins, 'non_sr': non_sr_proteins}

# ============================================================
# ANALYSIS 3: Stress Granule vs Ribosomal
# ============================================================
def analyze_compartments():
    print("\n" + "=" * 70)
    print("ANALYSIS 3: Stress Granule vs Ribosomal Proteins")
    print("=" * 70)

    stress_granule_keywords = ['stress granule', 'P-body', 'processing body',
                               'cytoplasmic granule', 'RNP granule']
    ribosome_keywords = ['ribosome', 'ribosomal', 'translation']

    all_accessions = list(set(d['accession'] for d in domains))[:100]

    sg_proteins = []
    ribosomal_proteins = []
    other_proteins = []

    print(f"\nClassifying {len(all_accessions)} proteins by compartment...")

    for i, acc in enumerate(all_accessions):
        info = get_protein_go_and_features(acc)
        combined = ' '.join(info['go_terms'] + info['locations'] + [info['name']]).lower()

        # get D/E
        acc_doms = [d for d in domains if d['accession'] == acc]
        de_pcts = [calculate_de(d['sequence']) for d in acc_doms]
        mean_de = np.mean(de_pcts) if de_pcts else 0

        is_sg = any(kw.lower() in combined for kw in stress_granule_keywords)
        is_ribo = any(kw.lower() in combined for kw in ribosome_keywords)

        protein_data = {'acc': acc, 'de': mean_de, 'name': info['name'][:50]}

        if is_sg:
            sg_proteins.append(protein_data)
        elif is_ribo:
            ribosomal_proteins.append(protein_data)
        else:
            other_proteins.append(protein_data)

        if (i+1) % 25 == 0:
            print(f"  {i+1}/{len(all_accessions)}")
        time.sleep(0.1)

    print(f"\nResults:")
    print(f"  Stress granule proteins: {len(sg_proteins)}")
    print(f"  Ribosomal proteins: {len(ribosomal_proteins)}")
    print(f"  Other proteins: {len(other_proteins)}")

    if sg_proteins:
        print(f"\n  Stress granule D/E: {np.mean([p['de'] for p in sg_proteins]):.2f}%")
    if ribosomal_proteins:
        print(f"  Ribosomal D/E: {np.mean([p['de'] for p in ribosomal_proteins]):.2f}%")

    if sg_proteins and ribosomal_proteins and len(sg_proteins) >= 3 and len(ribosomal_proteins) >= 3:
        stat, pval = stats.mannwhitneyu([p['de'] for p in sg_proteins],
                                        [p['de'] for p in ribosomal_proteins])
        print(f"\n  SG vs Ribosomal: p = {pval:.4f}")
        if pval < 0.05:
            print("  *** STRESS GRANULE PROTEINS HAVE DIFFERENT D/E ***")

    return {'sg': sg_proteins, 'ribosomal': ribosomal_proteins, 'other': other_proteins}

# ============================================================
# ANALYSIS 4: KH Regulatory Mechanisms
# ============================================================
def analyze_kh_mechanisms():
    print("\n" + "=" * 70)
    print("ANALYSIS 4: KH Domain Regulatory Mechanisms")
    print("=" * 70)

    # KH domains don't use D/E autoinhibition
    # what DO they use?

    kh_accessions = list(set(d['accession'] for d in domains if d['pfam_name'] == 'KH_1'))[:30]

    print(f"\nAnalyzing {len(kh_accessions)} KH proteins for regulatory features...")

    kh_data = []

    for i, acc in enumerate(kh_accessions):
        info = get_protein_go_and_features(acc)

        # get sequence composition
        acc_doms = [d for d in domains if d['accession'] == acc and d['pfam_name'] == 'KH_1']
        all_seq = ''.join(d['sequence'] for d in acc_doms)

        # check for various regulatory motifs
        de_pct = calculate_de(all_seq)
        kr_pct = calculate_kr(all_seq)
        gly_pct = all_seq.count('G') / len(all_seq) * 100 if all_seq else 0
        pro_pct = all_seq.count('P') / len(all_seq) * 100 if all_seq else 0

        # RGG motifs
        import re
        rgg_count = len(re.findall(r'RG{1,2}', all_seq))

        kh_data.append({
            'acc': acc,
            'name': info['name'][:40],
            'de': de_pct,
            'kr': kr_pct,
            'gly': gly_pct,
            'pro': pro_pct,
            'rgg': rgg_count,
            'phospho': info['phospho_count'],
            'go': info['go_terms']
        })

        if (i+1) % 15 == 0:
            print(f"  {i+1}/{len(kh_accessions)}")
        time.sleep(0.1)

    print(f"\nKH protein characteristics:")
    print(f"  D/E: {np.mean([p['de'] for p in kh_data]):.2f}% (vs RRM ~14%)")
    print(f"  K/R: {np.mean([p['kr'] for p in kh_data]):.2f}%")
    print(f"  Gly: {np.mean([p['gly'] for p in kh_data]):.2f}%")
    print(f"  Pro: {np.mean([p['pro'] for p in kh_data]):.2f}%")
    print(f"  RGG motifs: {np.mean([p['rgg'] for p in kh_data]):.1f}")
    print(f"  Phospho sites: {np.mean([p['phospho'] for p in kh_data]):.1f}")

    # compare to RRM
    print("\n  Comparison to RRM proteins:")
    rrm_doms = [d for d in domains if d['pfam_name'] == 'RRM_1'][:100]
    rrm_gly = [d['sequence'].count('G') / len(d['sequence']) * 100 for d in rrm_doms]
    rrm_pro = [d['sequence'].count('P') / len(d['sequence']) * 100 for d in rrm_doms]

    kh_gly_vals = [p['gly'] for p in kh_data]
    kh_pro_vals = [p['pro'] for p in kh_data]

    stat, pval = stats.mannwhitneyu(kh_gly_vals, rrm_gly)
    print(f"  KH Gly ({np.mean(kh_gly_vals):.1f}%) vs RRM Gly ({np.mean(rrm_gly):.1f}%): p = {pval:.4f}")

    stat, pval = stats.mannwhitneyu(kh_pro_vals, rrm_pro)
    print(f"  KH Pro ({np.mean(kh_pro_vals):.1f}%) vs RRM Pro ({np.mean(rrm_pro):.1f}%): p = {pval:.4f}")

    # look at GO terms for regulatory hints
    all_go = []
    for p in kh_data:
        all_go.extend(p['go'])

    print("\n  Most common GO terms in KH proteins:")
    go_counter = Counter(all_go)
    for term, count in go_counter.most_common(10):
        print(f"    {count}: {term}")

    return kh_data

# ============================================================
# MAIN
# ============================================================
def main():
    print("=" * 70)
    print("DEEP FUNCTIONAL ANALYSIS")
    print("=" * 70)

    results = {}

    # run all analyses
    results['decay'] = analyze_mrna_decay()
    results['sr_phospho'] = analyze_sr_phosphorylation()
    results['compartments'] = analyze_compartments()
    results['kh'] = analyze_kh_mechanisms()

    # save results
    # convert to serializable format
    output = {
        'decay_analysis': {
            'n_decay': len(results['decay']['decay']),
            'n_splicing': len(results['decay']['splicing']),
            'decay_de_mean': np.mean([p['de'] for p in results['decay']['decay']]) if results['decay']['decay'] else None,
            'splicing_de_mean': np.mean([p['de'] for p in results['decay']['splicing']]) if results['decay']['splicing'] else None,
        },
        'sr_analysis': {
            'n_sr': len(results['sr_phospho']['sr']),
            'n_non_sr': len(results['sr_phospho']['non_sr']),
        },
        'compartment_analysis': {
            'n_sg': len(results['compartments']['sg']),
            'n_ribosomal': len(results['compartments']['ribosomal']),
        }
    }

    with open(RESULTS_DIR / "deep_functional_analysis.json", 'w') as f:
        json.dump(output, f, indent=2, default=str)

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    print("\nResults saved to deep_functional_analysis.json")

    return results

if __name__ == "__main__":
    main()
