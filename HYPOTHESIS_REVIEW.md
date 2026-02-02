# D/E Autoinhibition Hypothesis: Critical Review & Strengthened Analysis

**Date:** 2026-02-01
**Status:** REVISED - Key finding overturned

---

## Executive Summary

The original hypothesis that "RBDs have D/E-rich segments serving autoinhibitory functions" has been **partially invalidated** by proper methodological analysis:

| Original Claim | Revised Finding | Impact |
|----------------|-----------------|--------|
| RBDs have more D/E segments (51% vs 35%) | **WRONG** - consecutive D/E segments are LESS common in RBDs (6.7% vs 7.9%) | ⚠️ CRITICAL |
| D/E enrichment (+2.05%) | Still valid at domain level, but Simpson's paradox applies | ⚠️ MODERATE |
| PTM enrichment near D/E (p=0.034) | After correction: marginal (Bonferroni p=0.10) | ⚠️ WEAK |
| D/E removal pathogenic (1.35x) | **CONFIRMED** (p=0.007) | ✓ STRONG |
| C-terminal D/E bias (+2.26%) | **CONFIRMED** (p < 10^-148) | ✓ STRONG |

---

## The Critical Revision

### Old Method (Arbitrary)
```python
# sliding window 8aa, threshold 37.5%
segments = find_de_rich_segments(seq, window_size=8, threshold=0.375)
```
**Result:** RBDs 51.2% vs non-RBDs 34.5% → OR = 1.99

### New Method (Literature Standard)
```python
# ≥3 consecutive D/E residues (Wang et al. standard)
segments = find_consecutive_de_segments(seq, min_length=3)
```
**Result:** RBDs **6.7%** vs non-RBDs **7.9%** → OR = **0.83** (opposite direction!)

### Why This Matters

The old sliding window method **overcounted** by:
1. Allowing non-D/E residues within "segments"
2. Using an arbitrary 37.5% threshold
3. Counting overlapping windows

The correct interpretation: **RBDs have elevated D/E CONTENT (+2.05%) but FEWER consecutive D/E runs**. This suggests:

- Autoinhibition uses **dispersed acidic residues**, not consecutive runs
- The charge distribution matters more than clustering
- Known autoinhibitory cases (Hfq, FBF-2) also lack consecutive runs in our detection

---

## What IS Supported

### 1. D/E Removal is Pathogenic ✓

| Test | Result | p-value |
|------|--------|---------|
| D/E removal vs other | 13.2% vs 9.3% pathogenic | **0.007** |
| Amidation (E→Q, D→N) | 22.3% pathogenic | **3.2×10^-5** |
| Charge reversal | 14.4% pathogenic | 0.07 (NS) |

**Interpretation:** Removing negative charge from RBPs is functionally deleterious, consistent with autoinhibition hypothesis.

### 2. C-terminal D/E Enrichment ✓

| Position | RBD D+E | non-RBD D+E | Diff |
|----------|---------|-------------|------|
| C-terminal | 12.28% | 10.02% | **+2.26%** |

χ² = 673, p < 10^-148

**Interpretation:** Matches Hfq acidic tail model - C-terminal acidic regions serve autoinhibitory functions.

### 3. Flanking Regions Are Disordered ✓

| Region | pLDDT | Interpretation |
|--------|-------|----------------|
| Domain core | 90.1 | Ordered |
| N-terminal flank | 68.3 | **Disordered** |
| C-terminal flank | 69.8 | **Disordered** |

**Interpretation:** D/E autoinhibitory regions reside in disordered flanks, not within structured domains.

### 4. Conservation in Specific Families ✓

| Family | D/E Conservation vs non-D/E | D/E→D/E Maintained |
|--------|-----------------------------|--------------------|
| LSM | **+15.7%** | 67.5% |
| RRM_1 | **+8.1%** | 62.0% |
| DEAD | **+6.0%** | 56.7% |
| KH_1 | -8.9% | 45.9% |
| S1 | -8.2% | 33.8% |

**Interpretation:** D/E conservation is **family-specific**, not universal.

### 5. Mutation Predictions Consistent ✓

| Statistic | Value |
|-----------|-------|
| Mean predicted binding increase | 4.97x |
| Predicted pathogenicity | 18.7% (vs 9.8% overall) |
| Top family (WW) | 5.86x predicted increase |

**Interpretation:** Terminal D/E in high-density clusters have strongest predicted effects, consistent with autoinhibition model.

---

## What is NOT Supported

### 1. "RBDs Have More D/E Segments" ✗

- Old: 51.2% vs 34.5% (OR = 1.99)
- **New: 6.7% vs 7.9% (OR = 0.83)** ← OPPOSITE

### 2. Universal D/E Enrichment ✗

- Domain-level: p < 10^-213 ✓
- **Family-level: p = 0.234** ✗
- Simpson's paradox is real

### 3. PTM Enrichment (After Correction) ✗

- Raw p = 0.034 → Bonferroni p = 0.10 (marginal)
- BH-FDR: also NS

---

## Revised Hypothesis

**Original:**
> "RBDs contain D/E-rich IDRs that electrostatically mimic nucleic acids for autoinhibition"

**Revised:**
> "SOME RBD families (RRM_1, LSM, DEAD, PUF) use DISPERSED D/E residues (not consecutive runs) in flanking IDRs for autoinhibitory regulation. This is family-specific and position-dependent (C-terminal bias). D/E removal is pathogenic, supporting functional importance."

---

## Remaining Uncertainties

| Question | Status | Needed |
|----------|--------|--------|
| CodonFM signal | Underpowered (n=15) | Better CDS source |
| MD validation | Not done | Simulations |
| Direct functional test | Literature only | Wet lab mutagenesis |

---

## Files Generated

### Strengthening Analyses
- `results/strengthened_analysis.txt` - Fixed D/E detection
- `results/mutation_predictions.json` - D/E→A effect predictions
- `results/mutation_predictions.txt` - Summary report

### Previous Analyses (Still Valid)
- `results/clinvar_de_analysis.txt` - Pathogenicity analysis
- `results/ptm_de_analysis.txt` - PTM enrichment
- `results/conservation_de_analysis.txt` - Cross-species conservation

---

## Conclusion

The D/E autoinhibition hypothesis is **partially supported** but requires significant revision:

1. **SUPPORTED:** C-terminal enrichment, pathogenicity of D/E removal, family-specific conservation
2. **NOT SUPPORTED:** Consecutive D/E segment enrichment (actually depleted in RBDs)
3. **MECHANISM:** Dispersed charge distribution, not consecutive runs

The strongest evidence remains the **ClinVar pathogenicity data** showing D/E removal is 1.35x more pathogenic (p=0.007) and amidation is 2.3x more pathogenic (p=3.2×10^-5).

---

*Analysis completed: 2026-02-01*
*Methodological issues identified and addressed*
*Key finding revised based on proper literature-standard detection*
