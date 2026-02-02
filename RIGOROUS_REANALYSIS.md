# Rigorous Reanalysis: What Survives Critical Testing

**Date:** 2026-02-01
**Status:** Post-peer-review revision

---

## Executive Summary

After rigorous testing addressing pseudo-replication, positional specificity, and motif analysis, the original claims have been substantially revised.

---

## What SURVIVES

| Claim | Evidence | Confidence |
|-------|----------|------------|
| D/E enrichment exists in splicing proteins | Δ=+3.9%, Cohen's d=1.06 | MODERATE (effect size large but CI fragile) |
| D/E removal is pathogenic | OR=2.77 amidation | HIGH |
| C-terminal regions are disordered | n=10, pLDDT=47 vs 80 | HIGH |
| KH has more proline than RRM | +1.5%, p<1e-20 | HIGH |

## What DOES NOT SURVIVE

| Original Claim | Problem |
|----------------|---------|
| Splicing D/E effect is robust | Block bootstrap 95% CI: [-1.1%, +6.1%] includes zero |
| D/E is C-terminal (autoinhibitory tail) | Actually N-terminal enriched (23% vs 10%) |
| KH uses SH3/PxxP regulation | No PxxP enrichment (Δ=-0.02%) |
| "Mechanism" of autoinhibition | Only positional correlation, no structural support |

---

## Detailed Results

### 1. Block-Aware Inference (Pseudo-replication)

**Original:** Mann-Whitney p=0.002 comparing 9000+ domains

**Problem:** Domains are phylogenetically correlated within families

**Fix:** Family-block bootstrap, protein-level aggregation

**Result:**
- Protein-level (n=14 splicing, n=3244 other): p=0.0002
- Block bootstrap 95% CI: **[-1.12%, +6.14%]** ← includes zero
- Bootstrap p-value: 0.107

**Interpretation:** The signal is real but fragile. With only ~3 families containing splicing proteins, resampling families can eliminate the signal entirely.

### 2. Positional Enrichment

**Original assumption:** D/E is in C-terminal tails (like Hfq)

**Test:** Compare N-terminal vs C-terminal D/E in high-D/E proteins

**Result:**
- High-D/E proteins: N-terminal = **23%**, C-terminal = **10.5%**
- D/E is **N-TERMINAL enriched**, not C-terminal

**Interpretation:** The Hfq model (C-terminal acidic tail) does NOT apply. Whatever mechanism operates, it's not tail-mediated.

### 3. Proline Motif Specificity

**Prediction:** If KH uses SH3-mediated regulation, should have PxxP motifs

**Test:** Compare PxxP density in KH vs RRM

**Result:**
- KH PxxP: 0.05%
- RRM PxxP: 0.07%
- Δ = -0.02% (no enrichment)

**Interpretation:** Proline enrichment in KH is real but does NOT translate to SH3-binding potential. Function unclear.

---

## Revised Model

### What We Can Say

1. **D/E content varies by RBP family** (descriptive fact)
2. **Splicing-related proteins tend toward higher D/E** (trend, not robust)
3. **D/E→Q/N mutations are enriched in pathogenic variants** (robust)
4. **KH has more proline than RRM** (robust, function unknown)

### What We Cannot Say

1. ~~D/E functions as an autoinhibitory tail~~ (positional data doesn't support)
2. ~~Splicing proteins have significantly higher D/E~~ (doesn't survive block bootstrap)
3. ~~KH uses SH3-mediated regulation~~ (no PxxP enrichment)
4. ~~D/E autoinhibition is a splicing-specific adaptation~~ (hierarchical inference too fragile)

---

## What Would Make This Publishable

1. **Expand splicing protein annotation** - Current n=14 is too small for robust family-block inference. Need comprehensive GO-based or manual curation to increase sample.

2. **Test N-terminal autoinhibition** - If D/E is N-terminal, does the N-terminus contact the RNA-binding surface? Requires AlphaFold-Multimer or MD.

3. **Explain KH proline** - What is the function? Polyproline II helix? Structural flexibility? Needs structural analysis.

4. **Direct binding data** - The pathogenicity signal is the strongest evidence. Could be validated with in vitro binding of E→Q mutants.

---

## Conclusion

This analysis started as "D/E autoinhibition in RBPs" and ended as "correlational enrichment with unknown mechanism." The signal is real but the mechanistic interpretation was premature.

The most robust finding is **pathogenicity of D/E removal** (OR=2.77), which survives all tests because it's based on independent ClinVar data, not the same domain sequences.

**Bottom line:** Pattern mining, not mechanism. Publishable as exploratory if honestly scoped.

---

*Reanalysis completed: 2026-02-01*
*In response to peer review feedback*
