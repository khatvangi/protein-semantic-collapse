#!/usr/bin/env python3
"""
Robust statistical analysis addressing all weaknesses:
1. Confound control (domain length, family size)
2. Multiple testing correction (FDR)
3. Bootstrap confidence intervals
4. Organism bias check
5. Sensitivity analysis for arbitrary thresholds
"""

import json
import numpy as np
from collections import Counter, defaultdict
from pathlib import Path
from scipy import stats
from scipy.stats import mannwhitneyu, chi2_contingency, pearsonr, spearmanr
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings('ignore')

STANDARD_AAS = set("ACDEFGHIKLMNPQRSTVWY")
ACIDIC = set("DE")


def compute_de_fraction(seq):
    valid = [aa for aa in seq if aa in STANDARD_AAS]
    if not valid:
        return 0.0
    return sum(1 for aa in valid if aa in ACIDIC) / len(valid)


def has_de_segment(seq, window=8, threshold=0.375):
    if len(seq) < window:
        return False
    for i in range(len(seq) - window + 1):
        w = seq[i:i+window]
        de = sum(1 for aa in w if aa in ACIDIC)
        if de / window >= threshold:
            return True
    return False


def bootstrap_ci(data, stat_func=np.mean, n_boot=2000, ci=95):
    """Bootstrap confidence interval."""
    if len(data) < 2:
        return np.nan, np.nan, np.nan
    boot_stats = []
    data = np.array(data)
    for _ in range(n_boot):
        sample = np.random.choice(data, size=len(data), replace=True)
        boot_stats.append(stat_func(sample))
    lower = np.percentile(boot_stats, (100 - ci) / 2)
    upper = np.percentile(boot_stats, 100 - (100 - ci) / 2)
    return stat_func(data), lower, upper


def main():
    base_dir = Path("/storage/kiran-stuff/protein-semantic-collapse")

    print("Loading data...")
    rbd_domains = []
    nonrbd_domains = []
    all_domains = []

    with open(base_dir / "data/domain_sequences.jsonl") as f:
        for line in f:
            d = json.loads(line)
            d["de_frac"] = compute_de_fraction(d["sequence"])
            d["has_de_seg"] = has_de_segment(d["sequence"])
            d["length"] = len(d["sequence"])
            all_domains.append(d)
            if d["is_rbd"]:
                rbd_domains.append(d)
            else:
                nonrbd_domains.append(d)

    print(f"  Total: {len(all_domains)}")
    print(f"  RBD: {len(rbd_domains)}, non-RBD: {len(nonrbd_domains)}")

    # =========================================================
    # CHECK 1: ORGANISM BIAS
    # =========================================================
    print("\n" + "=" * 70)
    print("CHECK 1: ORGANISM BIAS")
    print("=" * 70)

    # load organism info from FASTA headers
    organism_counts = {"RBD": Counter(), "non-RBD": Counter()}

    for fasta_type, domains in [("RBD", rbd_domains), ("non-RBD", nonrbd_domains)]:
        # we don't have organism in domain data, so check accession patterns
        # UniProt accessions: human often P/Q, mouse often Q
        acc_prefixes = Counter(d["accession"][0] for d in domains)
        print(f"\n{fasta_type} accession prefixes (proxy for organism):")
        for prefix, count in acc_prefixes.most_common(5):
            print(f"  {prefix}: {count} ({count/len(domains)*100:.1f}%)")

    # =========================================================
    # CHECK 2: DOMAIN LENGTH CONFOUND
    # =========================================================
    print("\n" + "=" * 70)
    print("CHECK 2: DOMAIN LENGTH CONFOUND")
    print("=" * 70)

    rbd_lengths = [d["length"] for d in rbd_domains]
    nonrbd_lengths = [d["length"] for d in nonrbd_domains]

    stat, p_length = mannwhitneyu(rbd_lengths, nonrbd_lengths)
    print(f"\nDomain length comparison:")
    print(f"  RBD mean:     {np.mean(rbd_lengths):.1f} (median: {np.median(rbd_lengths):.0f})")
    print(f"  non-RBD mean: {np.mean(nonrbd_lengths):.1f} (median: {np.median(nonrbd_lengths):.0f})")
    print(f"  Mann-Whitney p = {p_length:.2e}")

    # correlation: length vs D/E
    all_lengths = [d["length"] for d in all_domains]
    all_de = [d["de_frac"] for d in all_domains]
    r_len_de, p_len_de = pearsonr(all_lengths, all_de)
    print(f"\nLength-D/E correlation: r = {r_len_de:.3f}, p = {p_len_de:.2e}")

    if abs(r_len_de) > 0.1 and p_len_de < 0.05:
        print("  ⚠️  LENGTH IS A CONFOUND - need to control for it")
    else:
        print("  ✓ Length is not strongly correlated with D/E")

    # =========================================================
    # CHECK 3: FAMILY SIZE IMBALANCE
    # =========================================================
    print("\n" + "=" * 70)
    print("CHECK 3: FAMILY SIZE IMBALANCE")
    print("=" * 70)

    family_sizes = Counter(d["pfam"] for d in rbd_domains)
    print(f"\nRBD family size distribution:")
    print(f"  Total families: {len(family_sizes)}")
    print(f"  Largest: {family_sizes.most_common(1)[0]}")
    print(f"  Smallest (>=10): {min(c for c in family_sizes.values() if c >= 10)}")
    print(f"  Top 5 families account for: {sum(c for _, c in family_sizes.most_common(5))/len(rbd_domains)*100:.1f}% of domains")

    # =========================================================
    # ROBUST TEST 1: D/E CONTENT WITH LENGTH CONTROL (Regression)
    # =========================================================
    print("\n" + "=" * 70)
    print("ROBUST TEST 1: D/E CONTENT (Regression with confound control)")
    print("=" * 70)

    # prepare data for regression
    y = np.array([d["de_frac"] for d in all_domains])
    X = np.column_stack([
        [1 if d["is_rbd"] else 0 for d in all_domains],  # is_rbd
        [d["length"] for d in all_domains],               # length
        [np.log(d["length"]) for d in all_domains],       # log_length
    ])
    X = sm.add_constant(X)

    model = sm.OLS(y, X).fit()

    print(f"\nOLS Regression: D/E ~ is_RBD + length + log_length")
    print(f"\n{'Variable':<15} {'Coef':>10} {'Std Err':>10} {'t':>8} {'p-value':>12}")
    print("-" * 60)
    var_names = ["intercept", "is_RBD", "length", "log_length"]
    for i, name in enumerate(var_names):
        print(f"{name:<15} {model.params[i]:>10.5f} {model.bse[i]:>10.5f} {model.tvalues[i]:>8.2f} {model.pvalues[i]:>12.2e}")

    print(f"\nR² = {model.rsquared:.4f}")
    print(f"\nis_RBD effect after controlling for length:")
    print(f"  Coefficient: {model.params[1]:.4f} ({model.params[1]*100:.2f}% D/E)")
    print(f"  p-value: {model.pvalues[1]:.2e}")

    if model.pvalues[1] < 0.001:
        print(f"  ✓ CONFIRMED: RBD effect is significant after length control")
    else:
        print(f"  ✗ NOT significant after controlling for length")

    # =========================================================
    # ROBUST TEST 2: D/E SEGMENT WITH LENGTH CONTROL (Logistic)
    # =========================================================
    print("\n" + "=" * 70)
    print("ROBUST TEST 2: D/E SEGMENT PRESENCE (Logistic regression)")
    print("=" * 70)

    y_seg = np.array([1 if d["has_de_seg"] else 0 for d in all_domains])

    try:
        logit_model = sm.Logit(y_seg, X).fit(disp=0)

        print(f"\nLogistic Regression: has_DE_segment ~ is_RBD + length + log_length")
        print(f"\n{'Variable':<15} {'Coef':>10} {'Std Err':>10} {'z':>8} {'p-value':>12} {'OR':>8}")
        print("-" * 70)
        for i, name in enumerate(var_names):
            odds_ratio = np.exp(logit_model.params[i]) if i > 0 else np.nan
            print(f"{name:<15} {logit_model.params[i]:>10.4f} {logit_model.bse[i]:>10.4f} {logit_model.tvalues[i]:>8.2f} {logit_model.pvalues[i]:>12.2e} {odds_ratio:>8.2f}")

        print(f"\nis_RBD effect (odds ratio): {np.exp(logit_model.params[1]):.2f}")
        print(f"p-value: {logit_model.pvalues[1]:.2e}")

        if logit_model.pvalues[1] < 0.001:
            print(f"  ✓ CONFIRMED: RBD effect on D/E segments significant after length control")
    except Exception as e:
        print(f"  Logistic regression failed: {e}")

    # =========================================================
    # ROBUST TEST 3: SENSITIVITY ANALYSIS (Different thresholds)
    # =========================================================
    print("\n" + "=" * 70)
    print("ROBUST TEST 3: SENSITIVITY ANALYSIS (Threshold variation)")
    print("=" * 70)

    print(f"\nD/E segment detection with different thresholds:")
    print(f"{'Window':<8} {'Thresh':<8} {'RBD %':>10} {'non-RBD %':>12} {'OR':>8} {'p-value':>12}")
    print("-" * 65)

    results_sensitivity = []
    for window in [6, 8, 10, 12]:
        for thresh in [0.3, 0.375, 0.4, 0.5]:
            rbd_pos = sum(1 for d in rbd_domains if has_de_segment(d["sequence"], window, thresh))
            nonrbd_pos = sum(1 for d in nonrbd_domains if has_de_segment(d["sequence"], window, thresh))

            rbd_pct = rbd_pos / len(rbd_domains) * 100
            nonrbd_pct = nonrbd_pos / len(nonrbd_domains) * 100

            table = [[rbd_pos, len(rbd_domains) - rbd_pos],
                     [nonrbd_pos, len(nonrbd_domains) - nonrbd_pos]]
            try:
                chi2, p_val, _, _ = chi2_contingency(table)
                or_val = (rbd_pos / (len(rbd_domains) - rbd_pos)) / (nonrbd_pos / (len(nonrbd_domains) - nonrbd_pos)) if nonrbd_pos > 0 else np.inf
            except:
                p_val = 1.0
                or_val = np.nan

            results_sensitivity.append((window, thresh, rbd_pct, nonrbd_pct, or_val, p_val))

            if window == 8:  # only print for window=8 to keep output manageable
                print(f"{window:<8} {thresh:<8.3f} {rbd_pct:>9.1f}% {nonrbd_pct:>11.1f}% {or_val:>8.2f} {p_val:>12.2e}")

    # check if all are significant
    all_sig = all(r[5] < 0.001 for r in results_sensitivity)
    print(f"\nAll threshold combinations significant (p<0.001): {'✓ YES' if all_sig else '✗ NO'}")

    # =========================================================
    # ROBUST TEST 4: BOOTSTRAP WHOLE ANALYSIS
    # =========================================================
    print("\n" + "=" * 70)
    print("ROBUST TEST 4: BOOTSTRAP ANALYSIS (n=1000)")
    print("=" * 70)

    n_boot = 1000
    boot_de_diffs = []
    boot_seg_ors = []

    for i in range(n_boot):
        # resample with replacement
        rbd_sample = [rbd_domains[j] for j in np.random.randint(0, len(rbd_domains), len(rbd_domains))]
        nonrbd_sample = [nonrbd_domains[j] for j in np.random.randint(0, len(nonrbd_domains), len(nonrbd_domains))]

        # D/E difference
        rbd_de = np.mean([d["de_frac"] for d in rbd_sample])
        nonrbd_de = np.mean([d["de_frac"] for d in nonrbd_sample])
        boot_de_diffs.append(rbd_de - nonrbd_de)

        # segment OR
        rbd_seg = sum(d["has_de_seg"] for d in rbd_sample)
        nonrbd_seg = sum(d["has_de_seg"] for d in nonrbd_sample)

        rbd_noseg = len(rbd_sample) - rbd_seg
        nonrbd_noseg = len(nonrbd_sample) - nonrbd_seg

        if rbd_noseg > 0 and nonrbd_seg > 0 and nonrbd_noseg > 0:
            or_val = (rbd_seg / rbd_noseg) / (nonrbd_seg / nonrbd_noseg)
            boot_seg_ors.append(or_val)

    # D/E difference CI
    de_diff_mean = np.mean(boot_de_diffs)
    de_diff_ci = (np.percentile(boot_de_diffs, 2.5), np.percentile(boot_de_diffs, 97.5))
    de_diff_excludes_zero = de_diff_ci[0] > 0 or de_diff_ci[1] < 0

    # OR CI
    or_mean = np.mean(boot_seg_ors)
    or_ci = (np.percentile(boot_seg_ors, 2.5), np.percentile(boot_seg_ors, 97.5))
    or_excludes_one = or_ci[0] > 1 or or_ci[1] < 1

    print(f"\nD/E content difference:")
    print(f"  Mean: {de_diff_mean*100:+.2f}%")
    print(f"  95% CI: [{de_diff_ci[0]*100:+.2f}%, {de_diff_ci[1]*100:+.2f}%]")
    print(f"  CI excludes zero: {'✓ YES' if de_diff_excludes_zero else '✗ NO'}")

    print(f"\nD/E segment odds ratio:")
    print(f"  Mean OR: {or_mean:.2f}")
    print(f"  95% CI: [{or_ci[0]:.2f}, {or_ci[1]:.2f}]")
    print(f"  CI excludes 1.0: {'✓ YES' if or_excludes_one else '✗ NO'}")

    # =========================================================
    # ROBUST TEST 5: FDR CORRECTION FOR ALL TESTS
    # =========================================================
    print("\n" + "=" * 70)
    print("ROBUST TEST 5: MULTIPLE TESTING CORRECTION (FDR)")
    print("=" * 70)

    # collect all p-values from various tests
    all_pvalues = {
        "D/E content (regression)": model.pvalues[1],
        "D/E segment (logistic)": logit_model.pvalues[1] if 'logit_model' in dir() else 1.0,
        "Length difference": p_length,
        "Length-D/E correlation": p_len_de,
    }

    # add position-specific tests
    for i, pos_name in enumerate(["N-term", "Early-mid", "Middle", "Late-mid", "C-term"]):
        # quick chi-square for each position
        def get_pos_de(domains, bin_idx, n_bins=5):
            de, total = 0, 0
            for d in domains:
                seq = d["sequence"]
                n = len(seq)
                if n < n_bins:
                    continue
                start = int(bin_idx / n_bins * n)
                end = int((bin_idx + 1) / n_bins * n)
                for aa in seq[start:end]:
                    if aa in STANDARD_AAS:
                        total += 1
                        if aa in ACIDIC:
                            de += 1
            return de, total - de

        rbd_de, rbd_other = get_pos_de(rbd_domains, i)
        nonrbd_de, nonrbd_other = get_pos_de(nonrbd_domains, i)
        table = [[rbd_de, rbd_other], [nonrbd_de, nonrbd_other]]
        try:
            _, p_pos, _, _ = chi2_contingency(table)
        except:
            p_pos = 1.0
        all_pvalues[f"Position {pos_name}"] = p_pos

    # FDR correction
    pval_list = list(all_pvalues.values())
    names_list = list(all_pvalues.keys())
    reject, pvals_corrected, _, _ = multipletests(pval_list, method='fdr_bh')

    print(f"\n{'Test':<30} {'Raw p':>12} {'FDR p':>12} {'Sig':>6}")
    print("-" * 65)
    for name, raw_p, fdr_p, sig in zip(names_list, pval_list, pvals_corrected, reject):
        sig_marker = "✓" if sig else ""
        print(f"{name:<30} {raw_p:>12.2e} {fdr_p:>12.2e} {sig_marker:>6}")

    # =========================================================
    # ROBUST TEST 6: FAMILY-WEIGHTED ANALYSIS
    # =========================================================
    print("\n" + "=" * 70)
    print("ROBUST TEST 6: FAMILY-WEIGHTED ANALYSIS")
    print("=" * 70)

    # compute per-family means, then compare
    rbd_family_de = defaultdict(list)
    for d in rbd_domains:
        rbd_family_de[d["pfam"]].append(d["de_frac"])

    nonrbd_family_de = defaultdict(list)
    for d in nonrbd_domains:
        nonrbd_family_de[d["pfam"]].append(d["de_frac"])

    # family-level means
    rbd_fam_means = [np.mean(v) for v in rbd_family_de.values() if len(v) >= 10]
    nonrbd_fam_means = [np.mean(v) for v in nonrbd_family_de.values() if len(v) >= 10]

    print(f"\nFamily-level analysis (families with n>=10):")
    print(f"  RBD families: {len(rbd_fam_means)}")
    print(f"  non-RBD families: {len(nonrbd_fam_means)}")

    if rbd_fam_means and nonrbd_fam_means:
        stat, p_fam = mannwhitneyu(rbd_fam_means, nonrbd_fam_means)
        print(f"\nFamily-level D/E comparison:")
        print(f"  RBD mean of family means: {np.mean(rbd_fam_means)*100:.2f}%")
        print(f"  non-RBD mean of family means: {np.mean(nonrbd_fam_means)*100:.2f}%")
        print(f"  Mann-Whitney p = {p_fam:.2e}")

        if p_fam < 0.05:
            print(f"  ✓ CONFIRMED at family level")
        else:
            print(f"  ⚠️ NOT significant at family level")

    # =========================================================
    # FINAL ROBUST SUMMARY
    # =========================================================
    print("\n" + "=" * 70)
    print("FINAL ROBUST SUMMARY")
    print("=" * 70)

    print(f"""
FINDING                          ROBUST TEST                    VERDICT
--------------------------------------------------------------------------------
D/E content difference           Regression (length-controlled) {'✓' if model.pvalues[1] < 0.001 else '✗'}
                                 Bootstrap CI excludes 0        {'✓' if de_diff_excludes_zero else '✗'}
                                 Family-weighted                {'✓' if p_fam < 0.05 else '✗' if 'p_fam' in dir() else '?'}

D/E segment enrichment           Logistic (length-controlled)   {'✓' if logit_model.pvalues[1] < 0.001 else '✗'}
                                 Bootstrap OR CI excludes 1     {'✓' if or_excludes_one else '✗'}
                                 Sensitivity (all thresholds)   {'✓' if all_sig else '✗'}

Multiple testing                 FDR correction applied         ✓
Confounds                        Length controlled              ✓
Robustness                       Bootstrap (n=1000)             ✓
""")

    # final verdict
    n_confirmed = sum([
        model.pvalues[1] < 0.001,
        de_diff_excludes_zero,
        logit_model.pvalues[1] < 0.001 if 'logit_model' in dir() else False,
        or_excludes_one,
        all_sig,
    ])

    print(f"ROBUSTNESS SCORE: {n_confirmed}/5 tests passed")
    if n_confirmed >= 4:
        print("→ FINDINGS ARE ROBUST")
    elif n_confirmed >= 3:
        print("→ FINDINGS ARE MODERATELY ROBUST")
    else:
        print("→ FINDINGS MAY NOT BE ROBUST")


if __name__ == "__main__":
    main()
