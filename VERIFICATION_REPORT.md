# Verification Report: D/E Autoinhibition Analysis

**Date:** 2026-02-01
**Purpose:** Critical verification of all claims before acceptance

---

## Claims Tested

### Claim 1: Non-D/E Families Use Glycine-Proline Mechanism

**Original:** Non-D/E families have +1.88% more glycine (p=1.58e-41)

**Verification:**

| Check | Result |
|-------|--------|
| zf-CCHC outlier (12.3% Gly) | Removing it: still significant |
| Length confound | **FOUND**: D/E families are longer (92.6 vs 51.4 aa) |
| Gly correlates with length | r=-0.273, p<1e-40 |
| Length-controlled (50-80aa) | Difference: +0.42% (p=0.0037) |
| Regression residuals | Difference: +0.68% (p=1.32e-06) |

**VERDICT:** ✓ PARTIALLY VERIFIED
- Effect is REAL but SMALLER than reported
- Original +1.88% inflated by length confound
- True effect: **+0.4% to +0.7%** after controlling for length

---

### Claim 2: D/E Autoinhibition is Eukaryotic Innovation

**Original:** LUCA families have higher D/E (12.10% vs 10.49%, p=0.088)

**Problems Found:**
1. Family origin classifications were WRONG
2. Literature verification showed DEAD, dsRBD, Helicase_C are LUCA-origin (not eukaryotic)

**Corrected Analysis:**

| Origin | n families | Mean D/E | D/E-enriched |
|--------|------------|----------|--------------|
| LUCA | 7 | 11.57% | 3 (LSM, DEAD, Helicase_C) |
| Eukaryotic | 5 | 10.33% | 1 (RRM_1) |
| p-value | | 0.43 (NS) | Fisher p=0.58 (NS) |

**VERDICT:** ✗ NOT VERIFIED
- D/E content NOT different between LUCA and eukaryotic families
- D/E autoinhibition exists in BOTH ancient and eukaryotic families
- NOT an evolutionary innovation

---

### Claim 3: D/E Content is "Tunable" (High CV)

**Original:** All families have CV > 0.2, suggesting D/E is specially variable

**Verification:**

| Family | D/E CV | K/R CV | Gly CV | Leu CV | D/E unusual? |
|--------|--------|--------|--------|--------|--------------|
| RRM_1 | 0.293 | 0.268 | 0.314 | 0.414 | No |
| KH_1 | 0.305 | 0.238 | 0.247 | 0.515 | No |
| LSM | 0.260 | 0.293 | 0.288 | 0.241 | No |
| DEAD | 0.162 | 0.158 | 0.256 | 0.212 | No |

**VERDICT:** ✗ NOT VERIFIED
- D/E CV is NOT higher than other amino acids
- Leucine often has HIGHER CV than D/E
- D/E variation is TYPICAL, not specially "tunable"

---

## Summary Table

| Original Claim | Status | Corrected Finding |
|----------------|--------|-------------------|
| Gly-Pro +1.88% in non-D/E families | **PARTIALLY VERIFIED** | Effect is +0.5%, not +1.9% |
| D/E autoinhibition is eukaryotic | **NOT VERIFIED** | Present in both LUCA and eukaryotic families |
| D/E is specially tunable | **NOT VERIFIED** | CV is typical for amino acids |

---

## What We CAN Claim (Verified)

1. **Glycine is slightly elevated in non-D/E families** (+0.5% after length control, p<0.01)
   - This is a MINOR difference, not a major alternative mechanism

2. **D/E autoinhibition is NOT correlated with evolutionary origin**
   - Ancient families: LSM, DEAD, Helicase_C use it
   - Ancient families: KH, dsRBD, S1 do NOT use it
   - Eukaryotic families: RRM uses it
   - Eukaryotic families: zinc fingers do NOT use it

3. **D/E variation is NORMAL**
   - Not evidence for special "tunability"
   - Similar to variation in other amino acids

---

## Lessons Learned

1. **Always control for confounds** (length, family size, etc.)
2. **Verify classifications from literature** before analysis
3. **Compare to baseline** (is D/E CV actually high compared to other AAs?)
4. **Small sample sizes** (n=7 vs n=5 families) have low statistical power

---

## Literature Sources for Family Origins

| Family | Origin | Source |
|--------|--------|--------|
| KH | LUCA | [Wikipedia](https://en.wikipedia.org/wiki/KH_domain) |
| LSM/Hfq | LUCA | [PMC3710371](https://pmc.ncbi.nlm.nih.gov/articles/PMC3710371/) |
| DEAD-box | LUCA | [PMC3093544](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3093544/) |
| dsRBD | LUCA | [Nature Reviews](https://www.nature.com/articles/nrm1528) |
| RRM | Primarily Eukaryotic | [Wikipedia](https://en.wikipedia.org/wiki/RNA_recognition_motif) |

---

*Verification completed: 2026-02-01*
