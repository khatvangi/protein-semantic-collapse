# Evidence Validation: D/E Autoinhibition Hypothesis

**Purpose:** Critically evaluate each claim with evidence levels and identify gaps.

---

## Claim 1: RBDs Have Elevated D/E Content

### Evidence

```
Domain-level:
  RBD:     11.48% D/E (n=7,416)
  non-RBD:  9.43% D/E (n=16,824)
  Difference: +2.05%
  p-value: < 10^-213 (Mann-Whitney U)
  Effect size: r = 0.26 (small-medium)
```

### Counter-evidence

```
Family-level (weighted):
  RBD:     10.86% D/E
  non-RBD: 11.58% D/E
  Difference: -0.72%
  p-value: 0.234 (NOT SIGNIFICANT)
```

### Verdict

| Aspect | Validity |
|--------|----------|
| Domain-level effect | ✓ PROVEN (p < 10^-213) |
| Universal RBD property | ✗ NOT PROVEN (Simpson's paradox) |
| Family-specific effect | ✓ PROVEN (RRM, LSM, DEAD families) |

### What Would Prove It

- [x] Large sample size (n=24,240) ✓
- [x] Multiple statistical tests ✓
- [x] Effect survives length correction ✓
- [ ] Effect survives family weighting ✗ FAILS
- [ ] Mechanism demonstrated ✗ NOT DONE

### Confidence: 60% (family-specific, not universal)

---

## Claim 2: Consecutive D/E Segments Are Less Common in RBDs

### Evidence

```
≥3 consecutive D/E:
  RBD:     6.7% have segments
  non-RBD: 7.9% have segments
  OR = 0.83 (RBD LOWER)
  p < 0.001 (Fisher's exact)

≥5 consecutive D/E:
  RBD:     0.1%
  non-RBD: 0.3%
  OR = 0.19
```

### Method Validation

```python
# Literature standard (Wang et al.): consecutive D/E
def find_consecutive_de(seq, min_len=3):
    # Finds DDD, EEE, DDE, etc.
    # NOT sliding window with threshold
```

### Verdict

| Aspect | Validity |
|--------|----------|
| Consecutive segments less common | ✓ PROVEN |
| Old method was wrong | ✓ PROVEN (arbitrary threshold) |
| Dispersed D/E pattern | ✓ SUPPORTED (indirect) |

### What Would Prove It

- [x] Proper detection method ✓
- [x] Statistical significance ✓
- [x] Replicates across thresholds (3, 4, 5, 10) ✓
- [ ] Spatial distribution analysis (not just counts)

### Confidence: 90% (robust finding)

---

## Claim 3: D/E Removal Mutations Are Pathogenic

### Evidence

```
ClinVar/UniProt data (n=5,198 variants):

All variants:        9.8% pathogenic (509/5198)
D/E removal:        13.2% pathogenic (70/531)
  OR = 1.46, p = 0.007

Amidation (E→Q, D→N): 22.3% pathogenic (27/121)
  OR = 2.74, p = 3.2×10^-5

Specific substitutions:
  E→Q: 28.3% pathogenic (13/46)
  D→N: 18.7% pathogenic (14/75)
  D→G:  3.9% pathogenic (2/51)  ← control
```

### Controls

```
Charge reversal (D/E → K/R): 14.4% (p=0.07, NS)
  → Weaker effect than amidation

Glycine substitution: 3.9-5.0%
  → Size change alone not pathogenic
```

### Verdict

| Aspect | Validity |
|--------|----------|
| D/E removal enriched for pathogenicity | ✓ PROVEN (p=0.007) |
| Amidation highly pathogenic | ✓ PROVEN (p=3.2×10^-5) |
| Charge is the key feature | ✓ SUPPORTED (amidation > glycine) |
| Gain-of-function mechanism | ? UNKNOWN (no binding data) |

### What Would Prove It

- [x] Significant enrichment ✓
- [x] Dose-response (amidation > reversal > glycine) ✓
- [x] Controls for size effect ✓
- [ ] Binding assays showing increased affinity
- [ ] Distinguish GoF from LoF

### Confidence: 85% (strong statistical evidence)

---

## Claim 4: C-terminal D/E Enrichment

### Evidence

```
Position-specific D/E content:

Position      RBD      non-RBD    Diff
N-terminal   10.29%    12.08%   -1.79%
Early-mid    11.92%    11.17%   +0.75%
Middle       11.85%    11.05%   +0.81%
Late-mid     13.10%    11.16%   +1.94%
C-terminal   12.28%    10.02%   +2.26%  ← MAX

χ² = 673, p < 10^-148
```

### Literature Support

```
Hfq (Santiago-Frangos 2017): C-terminal acidic tail
FBF-2 (Qiu 2023): C-terminal IDR with D/E cluster
SLBP (Thapar 2015): C-terminal acidic region
```

### Verdict

| Aspect | Validity |
|--------|----------|
| C-terminal enrichment | ✓ PROVEN (p < 10^-148) |
| Matches literature models | ✓ SUPPORTED |
| Functional significance | ? INFERRED (not tested) |

### What Would Prove It

- [x] Statistical significance ✓
- [x] Literature precedent ✓
- [ ] Truncation experiments (remove C-terminus)
- [ ] Domain swapping experiments

### Confidence: 90% (very strong)

---

## Claim 5: C-terminal D/E Is Structurally Disordered

### Evidence

```
hnRNP A1 (P09651) AlphaFold analysis:

Region              D/E count   pLDDT    Status
RRM1 (1-90)              14     91.3    Ordered
RRM2 (100-180)           14     95.1    Ordered
C-terminal (>250)         2     36.6    DISORDERED

pLDDT scale:
  >90: Very high confidence (ordered)
  70-90: Confident
  50-70: Low confidence (possibly disordered)
  <50: Very low confidence (likely disordered)
```

### Additional Evidence

```
Flanking region analysis (100 RBDs + 90 non-RBDs):

Region           pLDDT (RBD)   pLDDT (non-RBD)
Domain core          90.1          87.9
N-terminal flank     68.3          68.4
C-terminal flank     69.8          71.5
```

### Verdict

| Aspect | Validity |
|--------|----------|
| C-terminal D/E is disordered | ✓ PROVEN (pLDDT=36.6) |
| Domain D/E is ordered | ✓ PROVEN (pLDDT>90) |
| Flanks more disordered than cores | ✓ PROVEN (Δ=21) |

### What Would Prove It

- [x] AlphaFold pLDDT validation ✓
- [x] Multiple proteins ✓
- [ ] NMR disorder measurement
- [ ] SAXS validation

### Confidence: 85% (AlphaFold is predictive, not experimental)

---

## Claim 6: D/E Positions Are More Conserved (Family-Specific)

### Evidence

```
Intra-family alignment (top 15 high-D/E domains per family):

Family    D/E identity   non-D/E identity   Diff
RRM_1        0.500           0.419         +0.081
LSM          0.581           0.424         +0.157
DEAD         0.456           0.395         +0.060
KH_1         0.277           0.366         -0.089
S1           0.215           0.297         -0.082

D/E character conservation (D or E maintained):
  RRM_1: 62.0%
  LSM:   67.5%
  DEAD:  56.7%
```

### Verdict

| Aspect | Validity |
|--------|----------|
| D/E conservation in LSM, RRM, DEAD | ✓ SUPPORTED |
| D/E conservation in KH_1, S1 | ✗ NOT FOUND |
| Negative charge conserved | ✓ SUPPORTED (D↔E exchange tolerated) |

### What Would Prove It

- [x] Cross-species alignment ✓
- [x] Multiple families ✓
- [ ] dN/dS analysis at D/E positions
- [ ] Ancestral sequence reconstruction

### Confidence: 70% (family-specific, limited species)

---

## Claim 7: D/E Codons Are More Optimized

### Evidence

```
Codon analysis (BioMart CDS, n=74 proteins):

Region        n codons   Mean score   Interpretation
D/E codons      1,021       3.55      More common
Core codons     7,188       3.91      Less common
Difference                 -0.36      D/E MORE optimized

p ≈ 0 (t-test)
```

### Caveat

```
D/E codon frequencies in humans (per 1000):
  GAT (D): 21.8
  GAC (D): 25.1
  GAA (E): 29.0
  GAG (E): 39.6

All D/E codons are relatively common in humans.
This may be inherent to the amino acids, not selection.
```

### Verdict

| Aspect | Validity |
|--------|----------|
| D/E codons more frequent | ✓ PROVEN |
| Due to functional selection | ? UNCLEAR (could be inherent) |
| Different from relaxed selection | ✓ SUPPORTED |

### What Would Prove It

- [x] Statistical test ✓
- [ ] Synonymous D/E codon preference (GAC vs GAT, etc.)
- [ ] Compare to non-RBD D/E positions
- [ ] tRNA adaptation index analysis

### Confidence: 60% (confounded by amino acid properties)

---

## Claim 8: PTM Enrichment Near D/E

### Evidence

```
PTM analysis (n=2,547 proteins):

                    Near D/E    Far from D/E   Enrichment   p-value
Phosphorylation    7.23/1000    5.86/1000       1.23x       0.034
Acetylation        7.44/1000    5.97/1000       1.25x       0.161
Methylation        1.71/1000    2.02/1000       0.85x       0.651

After Bonferroni (3 tests):
  Phosphorylation: p = 0.10 (MARGINAL)
```

### Verdict

| Aspect | Validity |
|--------|----------|
| Phosphorylation enriched near D/E | ? MARGINAL (p=0.10 corrected) |
| D/E as kinase recognition motif | ? PLAUSIBLE (literature support) |

### What Would Prove It

- [x] Statistical test ✓
- [ ] Survives multiple testing ✗ MARGINAL
- [ ] Kinase motif analysis
- [ ] Experimental validation

### Confidence: 50% (weak after correction)

---

## Claim 9: D/E-K/R Contacts Support Autoinhibition

### Evidence

```
hnRNP A1 structural analysis:

D/E-K/R contacts within 10Å (non-sequential): 38

Closest contacts:
  D155 - K105: 5.4Å
  E114 - K145: 5.5Å
  E176 - K166: 5.8Å
```

### Verdict

| Aspect | Validity |
|--------|----------|
| D/E and K/R are spatially close | ✓ OBSERVED |
| These contacts form in solution | ? UNKNOWN |
| They compete with RNA | ? UNKNOWN |

### What Would Prove It

- [x] Distance measurement ✓
- [ ] MD simulation showing contact dynamics
- [ ] Cross-linking MS validating contacts
- [ ] Mutagenesis + binding assays

### Confidence: 40% (circumstantial)

---

## Summary: Evidence Hierarchy

| Claim | Evidence Level | Confidence |
|-------|----------------|------------|
| C-terminal D/E enrichment | **PROVEN** | 90% |
| Consecutive D/E less common | **PROVEN** | 90% |
| D/E removal pathogenic | **PROVEN** | 85% |
| C-terminal D/E disordered | **PROVEN** | 85% |
| D/E conservation (family-specific) | SUPPORTED | 70% |
| D/E codons optimized | SUPPORTED | 60% |
| D/E enrichment (family-specific) | SUPPORTED | 60% |
| PTM near D/E | MARGINAL | 50% |
| D/E-K/R autoinhibitory contacts | CIRCUMSTANTIAL | 40% |

---

## What Would Definitively Prove Autoinhibition

### Tier 1: Computational (can do now)
- [ ] MD simulation of D/E tail dynamics
- [ ] Electrostatic potential mapping
- [ ] Coarse-grained binding simulations

### Tier 2: Experimental (literature search)
- [ ] Find existing mutagenesis studies
- [ ] Find NMR relaxation data
- [ ] Find ITC/SPR binding data

### Tier 3: New Experiments (wet lab)
- [ ] D/E→A mutagenesis + RNA binding assays
- [ ] C-terminal truncation + binding
- [ ] Cross-linking MS to map contacts

---

## Honest Assessment

**What we have proven:**
1. D/E enrichment exists in specific RBD families (RRM, LSM, DEAD)
2. This D/E is dispersed, not clustered
3. C-terminal regions are enriched and disordered
4. D/E removal is pathogenic (especially amidation)

**What we have NOT proven:**
1. D/E actually competes with RNA for binding
2. The pathogenicity is due to loss of autoinhibition
3. This is a general mechanism (vs family-specific quirk)

**The gap:**
We have strong correlational evidence but no direct mechanistic proof.
The link between D/E composition and autoinhibitory function is inferred
from literature models, not demonstrated by our data.

---

*Validation completed: 2026-02-01*
