#!/usr/bin/env python3
"""
Phylogenetic Analysis of D/E Content Evolution within RBD Families

For each family:
1. Align representative sequences with MAFFT
2. Build ML tree with IQ-TREE
3. Map D/E content onto tree tips
4. Test for phylogenetic signal in D/E content

This tells us: Is D/E content conserved within families, or variable?
"""

import subprocess
import json
from pathlib import Path
from collections import defaultdict
import re

RESULTS_DIR = Path("/storage/kiran-stuff/protein-semantic-collapse/results")
PHYLO_DIR = RESULTS_DIR / "phylogenetics"

# direct binary paths (avoiding conda activate issues in subprocess)
MAFFT_BIN = "/storage/kiran-stuff/blast_env/bin/mafft"
IQTREE_BIN = "/storage/kiran-stuff/blast_env/bin/iqtree3"

def calculate_de_content(seq):
    """calculate D/E percentage"""
    seq_clean = seq.replace('-', '')  # remove gaps
    if not seq_clean:
        return 0
    return sum(1 for aa in seq_clean if aa in 'DE') / len(seq_clean) * 100

def run_mafft(input_fasta, output_aln):
    """align sequences with MAFFT"""
    cmd = f"{MAFFT_BIN} --auto --quiet {input_fasta} > {output_aln}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  MAFFT error: {result.stderr}")
        return False
    return True

def run_iqtree(alignment_file, output_prefix):
    """run IQ-TREE for ML phylogeny"""
    # use automatic model selection, fast bootstrap
    cmd = (f"{IQTREE_BIN} -s {alignment_file} "
           f"--prefix {output_prefix} "
           f"-m TEST "  # automatic model selection
           f"-B 1000 "  # ultrafast bootstrap
           f"--quiet "
           f"-T AUTO "  # auto threads (changed from -nt)
           f"--redo")  # overwrite if exists

    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  IQ-TREE error: {result.stderr[:500]}")
        return False
    return True

def parse_fasta(fasta_file):
    """parse FASTA file, return dict of id -> sequence"""
    sequences = {}
    current_id = None
    current_seq = []

    with open(fasta_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_id:
            sequences[current_id] = ''.join(current_seq)

    return sequences

def extract_de_from_header(header):
    """extract D/E content from sequence header if encoded"""
    match = re.search(r'DE([\d.]+)', header)
    if match:
        return float(match.group(1))
    return None

def analyze_de_variation(fasta_file):
    """analyze D/E content variation within family"""
    seqs = parse_fasta(fasta_file)

    de_contents = []
    for header, seq in seqs.items():
        de = calculate_de_content(seq)
        de_contents.append(de)

    if not de_contents:
        return None

    import numpy as np
    return {
        'n_seqs': len(de_contents),
        'mean_de': np.mean(de_contents),
        'std_de': np.std(de_contents),
        'min_de': np.min(de_contents),
        'max_de': np.max(de_contents),
        'range_de': np.max(de_contents) - np.min(de_contents),
        'cv_de': np.std(de_contents) / np.mean(de_contents) if np.mean(de_contents) > 0 else 0
    }

def create_de_trait_file(fasta_file, output_file):
    """create trait file for ancestral state reconstruction"""
    seqs = parse_fasta(fasta_file)

    with open(output_file, 'w') as f:
        f.write("taxon\tDE_content\n")
        for header, seq in seqs.items():
            de = calculate_de_content(seq)
            # clean header for tree compatibility
            clean_header = header.replace(':', '_').replace(' ', '_')
            f.write(f"{clean_header}\t{de:.2f}\n")

def main():
    print("=" * 70)
    print("Phylogenetic Analysis of D/E Content Evolution")
    print("=" * 70)

    # find all FASTA files
    fasta_files = list(PHYLO_DIR.glob("*_representatives.fasta"))
    print(f"\nFound {len(fasta_files)} family FASTA files")

    results = {}

    for fasta_file in sorted(fasta_files):
        family = fasta_file.stem.replace('_representatives', '')
        print(f"\n{'='*60}")
        print(f"Processing {family}")
        print(f"{'='*60}")

        # analyze D/E variation before alignment
        variation = analyze_de_variation(fasta_file)
        if variation:
            print(f"  D/E variation: {variation['mean_de']:.1f}% ± {variation['std_de']:.1f}")
            print(f"  Range: {variation['min_de']:.1f}% - {variation['max_de']:.1f}%")
            print(f"  CV: {variation['cv_de']:.3f}")

        # step 1: align with MAFFT
        aln_file = PHYLO_DIR / f"{family}.aln"
        print(f"  Aligning with MAFFT...")
        if not run_mafft(fasta_file, aln_file):
            print(f"  FAILED: Could not align {family}")
            continue

        # step 2: build tree with IQ-TREE
        tree_prefix = PHYLO_DIR / family
        print(f"  Building ML tree with IQ-TREE...")
        if not run_iqtree(aln_file, tree_prefix):
            print(f"  FAILED: Could not build tree for {family}")
            continue

        # step 3: create trait file for D/E content
        trait_file = PHYLO_DIR / f"{family}_de_traits.tsv"
        create_de_trait_file(fasta_file, trait_file)
        print(f"  Created D/E trait file: {trait_file.name}")

        # check if tree was created
        tree_file = PHYLO_DIR / f"{family}.treefile"
        if tree_file.exists():
            print(f"  ✓ Tree created: {tree_file.name}")

            # read tree to get basic stats
            with open(tree_file) as f:
                tree_str = f.read().strip()
            print(f"  Tree length: {len(tree_str)} chars")
        else:
            print(f"  ✗ No tree file found")

        results[family] = {
            'variation': variation,
            'alignment': str(aln_file),
            'tree': str(tree_file) if tree_file.exists() else None,
            'traits': str(trait_file)
        }

    # summary
    print("\n" + "=" * 70)
    print("SUMMARY: D/E Variation Within Families")
    print("=" * 70)

    print(f"\n{'Family':<15} {'Mean D/E':>10} {'Std':>8} {'Range':>10} {'CV':>8}")
    print("-" * 55)

    for family, data in sorted(results.items()):
        if data['variation']:
            v = data['variation']
            print(f"{family:<15} {v['mean_de']:>10.1f}% {v['std_de']:>8.1f} {v['range_de']:>10.1f} {v['cv_de']:>8.3f}")

    # interpretation
    print("\n" + "=" * 70)
    print("INTERPRETATION")
    print("=" * 70)

    # find high-variation vs low-variation families
    cv_values = [(fam, data['variation']['cv_de'])
                 for fam, data in results.items()
                 if data['variation']]

    cv_values.sort(key=lambda x: x[1], reverse=True)

    print("\nHigh D/E variation (CV > 0.2) suggests:")
    print("  - D/E content is NOT under strong purifying selection")
    print("  - May be evolving neutrally or adaptively")

    print("\nLow D/E variation (CV < 0.1) suggests:")
    print("  - D/E content is conserved")
    print("  - Functional constraint on D/E composition")

    print("\nBy family:")
    for fam, cv in cv_values:
        if cv > 0.2:
            print(f"  {fam}: CV={cv:.3f} - HIGH variation")
        elif cv < 0.1:
            print(f"  {fam}: CV={cv:.3f} - LOW variation (conserved)")
        else:
            print(f"  {fam}: CV={cv:.3f} - moderate")

    # save results
    with open(RESULTS_DIR / "phylogenetic_analysis_summary.json", 'w') as f:
        json.dump(results, f, indent=2, default=str)

    print(f"\nResults saved to {RESULTS_DIR / 'phylogenetic_analysis_summary.json'}")
    print(f"Tree files in {PHYLO_DIR}/")

if __name__ == "__main__":
    main()
