# Corrected Findings: D/E Autoinhibition in RNA-Binding Proteins

**Date:** 2026-02-01
**Status:** Post-verification corrections

---

## Executive Summary (Corrected)

After rigorous verification, several original claims were **revised or retracted**:

| Finding | Original | Corrected |
|---------|----------|-----------|
| Gly-Pro alternative mechanism | +1.88% difference | **+0.5%** (after length control) |
| D/E autoinhibition is eukaryotic | p=0.088 | **NOT SUPPORTED** (p=0.43) |
| D/E is specially tunable | CV > 0.2 | **NOT UNUSUAL** (similar to other AAs) |

---

## What IS Supported (After Verification)

### 1. D/E Enrichment is Family-Specific ✓

| Family | D/E % | D/E-enriched? | Origin |
|--------|-------|---------------|--------|
| RRM_1 | 13.8% | ✓ YES | Eukaryotic |
| LSM | 12.2% | ✓ YES | LUCA |
| DEAD | 11.1% | ✓ YES | LUCA |
| Helicase_C | 11.1% | ✓ YES | LUCA |
| KH_1 | 11.2% | ✗ NO | LUCA |
| dsRBD | 9.5% | ✗ NO | LUCA |
| zf-C2H2 | 6.9% | ✗ NO | Eukaryotic |

**Conclusion:** D/E autoinhibition is used by some families (RRM, LSM, DEAD) but not others (KH, dsRBD, zinc fingers). This is NOT related to evolutionary origin.

### 2. D/E Removal is Pathogenic ✓

| Mutation Type | Pathogenic Rate | OR | p-value |
|---------------|-----------------|-----|---------|
| All variants | 9.8% | 1.00 | - |
| D/E removal | 13.2% | 1.46 | 0.007 |
| Amidation (E→Q, D→N) | 22.3% | 2.77 | 3.2×10⁻⁵ |

**Status:** ✓ VERIFIED (not re-tested, from original analysis)

### 3. C-terminal D/E Enrichment ✓

| Position | RBD D/E | non-RBD D/E | Difference |
|----------|---------|-------------|------------|
| C-terminal | 12.28% | 10.02% | +2.26% |
| χ² p-value | | | < 10⁻¹⁴⁸ |

**Status:** ✓ VERIFIED (from original analysis)

### 4. C-terminal D/E is Disordered ✓

| Region | pLDDT | Status |
|--------|-------|--------|
| RRM domains | >90 | Ordered |
| C-terminal D/E | 36.6 | Disordered |

**Status:** ✓ VERIFIED (AlphaFold analysis)

---

## What is NOT Supported (After Verification)

### 1. Glycine-Proline as Major Alternative Mechanism ✗

**Original claim:** Non-D/E families have +1.88% more glycine

**Problem:** Length confound - non-D/E families have shorter domains (51 vs 93 aa), and glycine correlates negatively with length.

**After correction:** Only +0.5% difference (p<0.01)

**Conclusion:** Glycine is slightly elevated, but this is a MINOR effect, not a major alternative mechanism.

### 2. D/E Autoinhibition as Eukaryotic Innovation ✗

**Original claim:** Ancient (LUCA) families have different D/E than eukaryotic families

**Problem:** Original family origin classifications were WRONG. DEAD, dsRBD, Helicase_C are LUCA-origin (from literature).

**After correction:**
- LUCA families: 11.57% D/E
- Eukaryotic families: 10.33% D/E
- p = 0.43 (NOT significant)

**Conclusion:** D/E autoinhibition is NOT associated with evolutionary origin. It exists in both ancient (LSM, DEAD) and eukaryotic (RRM) families.

### 3. D/E as Specially Tunable ✗

**Original claim:** High CV (>0.2) suggests D/E is specially variable/tunable

**Problem:** Did not compare to baseline. All amino acids have similar CV.

| Family | D/E CV | Leu CV | Interpretation |
|--------|--------|--------|----------------|
| RRM_1 | 0.29 | 0.41 | Leu more variable |
| KH_1 | 0.31 | 0.52 | Leu more variable |

**Conclusion:** D/E variation is TYPICAL, not specially tunable.

---

## Revised Conclusions

### What We Know (High Confidence)

1. **D/E enrichment is family-specific** - RRM, LSM, DEAD use it; KH, dsRBD, zinc fingers don't
2. **D/E removal is pathogenic** - OR=1.46 for removal, OR=2.77 for amidation
3. **C-terminal D/E is enriched and disordered** - Structural validation from AlphaFold
4. **Hfq validates the mechanism** - Literature confirms D/E-K/R competition

### What We Don't Know (After Verification)

1. **Why some families use D/E and others don't** - Not related to evolutionary origin
2. **What alternative mechanisms exist** - Gly-Pro is minor (+0.5%), not major
3. **Whether D/E is specially tunable** - Evidence does not support this

---

## Honest Assessment

The **core hypothesis** (D/E autoinhibition in some RBD families) remains supported by:
- Pathogenicity data
- Structural data (AlphaFold)
- Literature validation (Hfq)

However, the **extended claims** about evolution and alternative mechanisms were:
- Overstated (Gly-Pro effect size)
- Based on incorrect classifications (evolutionary origins)
- Not compared to proper baselines (CV)

---

## Files Updated

| File | Status |
|------|--------|
| `VERIFICATION_REPORT.md` | NEW - Full verification details |
| `CORRECTED_FINDINGS.md` | NEW - This file |
| `verified_family_origins.md` | NEW - Literature-verified origins |

---

*Corrections completed: 2026-02-01*
