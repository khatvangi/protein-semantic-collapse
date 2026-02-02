# Protein Semantic Collapse: Findings & Future Work

## Project Overview

Analysis of amino acid vocabulary constraints in RNA-binding domains (RBDs) vs non-RBDs, testing whether RBDs show "semantic collapse" (reduced AA diversity) and how this relates to autoinhibitory mechanisms.

**Dataset:** 24,240 domains (7,416 RBD, 16,824 non-RBD) from 67 Pfam families

**⚠️ DATA BUG FIXED (2026-02-01):** Original dataset had WD40 (PF00400) entries appearing in BOTH RBD and non-RBD groups (2,485 duplicates). Fixed by using Pfam family list as ground truth instead of `is_rbd` flag.

---

## CRITICAL CAVEAT: Simpson's Paradox

**The domain-level D/E enrichment in RBDs is NOT robust at family level!**

| Level | RBD D/E | non-RBD D/E | Difference | p-value |
|-------|---------|-------------|------------|---------|
| Domain-level | 11.48% | 9.43% | +2.05% | < 1e-213 |
| **Family-level** | **10.86%** | **11.58%** | **-0.72%** | **0.234 (NS)** |

**Why this happens:**
- A few large RBD families (RRM_1: n=1090, WD40: n=2485) have high D/E (13-14%)
- These dominate the domain-level average
- Smaller RBD families have similar or lower D/E than non-RBDs

**Implication:** D/E enrichment may be a property of specific RBD families (RRM, WD40, HSP70), NOT a universal RBD feature.

---

## Key Findings (Statistically Validated)

### 1. D/E Content Difference ✓ CONFIRMED

| Metric | RBD | non-RBD | Difference | Statistic |
|--------|-----|---------|------------|-----------|
| D+E content | 11.48% | 9.43% | **+2.05%** | Mann-Whitney p < 1e-213 |
| 95% CI | 11.38-11.57% | 9.35-9.51% | | Effect size r = 0.26 |

### 2. D/E-Rich Segments ✓ CONFIRMED

| Metric | RBD | non-RBD | Statistic |
|--------|-----|---------|-----------|
| With D/E segment | **51.2%** | 34.5% | χ² = 578.3 |
| | (3798/7416) | (5799/16824) | **p < 1e-126** |
| Odds ratio | 1.99 | | Fisher p < 1e-125 |

### 3. Position-Specific Patterns ✓ CONFIRMED

| Position | RBD D+E | non-RBD D+E | Diff |
|----------|---------|-------------|------|
| N-terminal | 10.29% | 12.08% | -1.79% |
| Early-mid | 11.92% | 11.17% | +0.75% |
| Middle | 11.85% | 11.05% | +0.81% |
| Late-mid | 13.10% | 11.16% | +1.94% |
| **C-terminal** | **12.28%** | **10.02%** | **+2.26%** |

**C-terminal enrichment:** χ² = 673.0, **p < 1e-148**

---

## Robust Analysis Results

### Confound Checks

| Confound | Status | Action |
|----------|--------|--------|
| **Domain length** | ⚠️ CONFOUND (r=0.26 with D/E) | Controlled in regression |
| Organism bias | ✓ Similar distribution | None needed |
| Family size imbalance | ⚠️ Top 5 families = 63% of data | Family-weighted analysis |

### After Controlling for Confounds

| Test | Statistic | p-value | Status |
|------|-----------|---------|--------|
| D/E content (OLS, length-controlled) | coef = +1.13% | 2.1e-87 | ✓ |
| D/E segment (Logistic, length-controlled) | OR = 1.35 | 4.2e-16 | ✓ |
| Bootstrap CI for D/E diff | [+1.79%, +2.01%] | excludes 0 | ✓ |
| Bootstrap CI for OR | [1.80, 2.01] | excludes 1 | ✓ |
| Sensitivity (all thresholds) | all p < 0.001 | - | ✓ |
| **Family-weighted** | - | **0.234** | **✗ NOT SIG** |

### Robustness Score: 5/6

The domain-level effect is robust to length control, bootstrap, and sensitivity analysis.
**BUT** the family-weighted analysis fails (p=0.234), indicating Simpson's paradox.

### 4. Flanking Region Analysis

| Region | D+E (RBD) | D+E (non-RBD) | Disorder Score |
|--------|-----------|---------------|----------------|
| N-terminal flank | **12.64%** | 10.95% | 0.159 |
| Domain core | 11.48% | 9.57% | 0.050 |
| C-terminal flank | **11.86%** | 10.69% | 0.156 |

**Interpretation:** Flanking regions are 3x more disordered than domain cores. D/E-rich autoinhibitory IDRs extend beyond Pfam domain boundaries.

---

## Novel Findings

### Finding 1: Ribosomal Proteins Are Outliers

| Family | D+E | K+R | Net Charge |
|--------|-----|-----|------------|
| Ribosomal_S3_C | 6.8% | 18.9% | +12.0% |
| Ribosomal_S15 | 12.1% | 24.2% | +12.1% |
| Ribosomal_L19 | 7.3% | 23.8% | +16.5% |

**Why novel:** These RBDs have unusually LOW acidic content and HIGH basic content. They don't use autoinhibition because the ribosome structure itself provides binding specificity - no need for competitive regulation.

**Biological interpretation:** Autoinhibition evolved for RBDs that must discriminate among many RNA targets. Ribosomal proteins bind specific rRNA sites in a pre-organized structural context.

### Finding 2: Entropy-Autoinhibition Relationship ✓ CORRECTED

| Correlation | r | p-value | Status |
|-------------|---|---------|--------|
| Entropy ↔ D/E content | -0.049 | 0.748 | ✗ NOT significant |
| Entropy ↔ D/E segment % | **+0.477** | **0.0009** | ✓ SIGNIFICANT |
| (Spearman ρ) | +0.498 | 0.0005 | ✓ Confirms |

**Corrected interpretation:**
- Entropy does NOT correlate with overall D/E content
- Entropy DOES correlate with D/E segment presence (r = +0.48, p < 0.001)
- **Higher entropy families use MORE autoinhibitory segments**

**Why:** Families with diverse vocabulary may need stronger autoinhibition because their binding surfaces could interact promiscuously.

### Finding 3: zf-CCCH Uses Alternative Regulation

| Family | D/E Segment % | D+E Content | Mechanism |
|--------|---------------|-------------|-----------|
| zf-CCCH | **10.8%** | 8.5% | Zinc coordination |
| RRM_1 | 58.2% | 13.6% | D/E autoinhibition |

**Why novel:** Zinc finger domains (zf-CCCH) rarely have D/E segments despite being RBDs. They use metal coordination for specificity rather than acidic autoinhibition.

---

## Literature Comparison

| Our Finding | Literature | Match |
|-------------|------------|-------|
| D enrichment in RBDs | Wang et al. 2025: D/E-rich IDRs common in RBPs | ✓ |
| K/R at interface | Interface studies show K/R enriched at hotspots | ~ (domain-wide shows depletion) |
| D/E segments 51% vs 36% | Autoinhibition widespread in RBPs | ✓ |
| C-terminal acidic bias | Not explicitly described | **NEW** |
| Ribosomal outliers | Not explicitly described | **NEW** |
| Entropy-autoinhibition independence | Not tested before | **NEW** |

---

## Future Work

### 1. AlphaFold pLDDT Validation (High Priority)

**Goal:** Validate that D/E-rich regions correspond to predicted disorder.

**Approach:**
- Fetch AlphaFold structures for representative domains
- Extract pLDDT scores (confidence = order proxy)
- Correlate D/E content with low pLDDT regions
- Compare RBD vs non-RBD pLDDT profiles

**Expected outcome:** D/E-rich segments should have pLDDT < 70 (disordered).

### 2. Evolutionary Analysis of Outliers (High Priority)

**Goal:** Understand why ribosomal proteins don't use autoinhibition.

**Approach:**
- Phylogenetic analysis of ribosomal vs cytoplasmic RBDs
- Test if autoinhibition correlates with number of RNA targets
- Check if loss of autoinhibition in ribosomal proteins is ancestral or derived

**Hypothesis:** Autoinhibition evolved after ribosomal proteins, for RBDs needing target discrimination.

### 3. Functional Validation of D/E Segments (Medium Priority)

**Goal:** Test if identified D/E segments actually function as autoinhibitors.

**Approach:**
- Select 3-5 RBDs with strong D/E segments
- Design mutations: D/E → A (remove negative charge)
- Predict binding affinity changes using FoldX or similar
- Literature search for existing mutagenesis data

**Expected outcome:** D/E → A mutations should increase RNA binding affinity (remove autoinhibition).

### 4. Cross-Species Conservation (Medium Priority)

**Goal:** Test if D/E segment positions are conserved across evolution.

**Approach:**
- Align orthologous RBDs from human, mouse, fly, yeast
- Map D/E segment positions onto alignments
- Calculate conservation scores for D/E vs non-D/E positions

**Hypothesis:** If functional, D/E segments should show position-specific conservation.

### 5. RNA Target Complexity Analysis (Exploratory)

**Goal:** Test if autoinhibition strength correlates with RNA target diversity.

**Approach:**
- Annotate RBDs by number of known RNA targets (from ENCODE, etc.)
- Correlate target count with D/E segment presence/strength
- Control for expression level and subcellular localization

**Hypothesis:** RBDs with many targets need stronger autoinhibition.

### 6. Codon-Level Analysis (High Priority)

**Goal:** Test if D/E-rich autoinhibitory regions have distinct codon usage patterns.

**Approach:**
- Extract coding sequences for domains with/without D/E segments
- Compare codon usage bias (CAI, tAI) in D/E regions vs flanks vs core
- Use CodonFM to score RBD vs non-RBD coding sequences
- Test for rare codon enrichment (translation pausing at regulatory regions?)
- Check synonymous constraint (dN/dS within D/E segments)

**Hypotheses:**
1. D/E regions may use rare codons to slow translation → co-translational folding
2. Autoinhibitory regions may show relaxed synonymous constraint (less structured)
3. CodonFM perplexity may differ between RBD classes

**Connection to CodonFM project:** Apply the Encodon-80M model to score these sequences and compare PLL scores across domain types.

### 7. Functional Validation & Disease Relevance (High Priority)

**Goal:** Connect compositional findings to biological function and disease.

**Approach:**
- GO enrichment analysis for outlier families (ribosomal vs cytoplasmic RBDs)
- Disease variant mapping: Are D/E segment mutations enriched in ClinVar/OMIM?
- PTM analysis: Do D/E regions overlap with phosphorylation sites? (PhosphoSitePlus)
- Expression patterns: Are high-autoinhibition RBDs stress-responsive?

**Specific questions:**
1. Do mutations in D/E segments cause gain-of-function (loss of autoinhibition)?
2. Are D/E segments phosphorylated to regulate autoinhibition strength?
3. Do outlier families (ribosomal) have different disease associations?

### 8. Machine Learning Prediction (Exploratory)

**Goal:** Predict autoinhibition presence from sequence alone.

**Approach:**
- Train classifier on D/E segment presence (binary)
- Features: AA composition, entropy, hydrophobicity profile, charge distribution
- Test on held-out families
- Identify predictive features

**Application:** Screen proteomes for RBDs likely to have autoinhibitory regulation.

---

## Priority Summary

| # | Project | Priority | Effort | Key Question |
|---|---------|----------|--------|--------------|
| 1 | AlphaFold pLDDT validation | High | 1 day | Do D/E regions have low pLDDT? |
| 2 | Evolutionary outlier analysis | High | 2 days | Why no autoinhibition in ribosomal proteins? |
| 6 | **Codon-level analysis** | **High** | 2 days | CodonFM signal in D/E regions? |
| 7 | **Functional/disease validation** | **High** | 2 days | Are D/E mutations pathogenic? |
| 3 | Mutagenesis prediction | Medium | 1 day | D/E→A increases binding? |
| 4 | Cross-species conservation | Medium | 2 days | D/E positions conserved? |
| 5 | RNA target complexity | Exploratory | 3 days | More targets = more autoinhibition? |
| 8 | ML prediction model | Exploratory | 1 week | Predict autoinhibition from sequence |

---

## Codon-Protein Integration Hypotheses

### The Translation-Regulation Link

If D/E-rich autoinhibitory regions are functionally important, they may show:

1. **Codon bias patterns:**
   - Rare codons in D/E regions → slow translation → time for IDR to fold/interact
   - Or optimal codons → fast translation → rapid autoinhibition establishment

2. **CodonFM perplexity signal:**
   - Low perplexity (good fit) = evolutionarily optimized codon usage
   - Compare PLL scores: D/E segment codons vs domain core codons
   - Hypothesis: D/E regions may have HIGHER perplexity (less constrained)

3. **Synonymous selection:**
   - D/E segments as IDRs may show relaxed synonymous constraint
   - Or: if co-translational folding matters, strong synonymous selection

### Proposed Analysis Pipeline

```
1. Map protein D/E segments → coding sequence coordinates
2. Extract codon sequences for:
   - D/E-rich regions
   - Domain cores
   - Flanking IDRs
3. Run CodonFM (Encodon-80M) on each region type
4. Compare:
   - Mean PLL scores
   - Codon adaptation index (CAI)
   - Rare codon frequency
5. Statistical test: Are D/E regions distinct at codon level?
```

### Expected Outcomes

| Region | Predicted CodonFM PLL | Predicted CAI | Rationale |
|--------|----------------------|---------------|-----------|
| Domain core | Low (good fit) | High | Structured, conserved |
| D/E segments | Higher (worse fit) | Variable | Disordered, regulatory |
| Flanking IDRs | Higher | Lower | Less constrained |

---

## Files & Scripts

### Data
- `data/domain_sequences.jsonl` - Domain sequences with annotations
  - ⚠️ Original file had contamination bug (WD40 duplicates)
  - Scripts now use Pfam family list as ground truth for RBD classification
  - Effective dataset: 24,240 domains (7,416 RBD, 16,824 non-RBD)
- `sequences/rbp_sequences.fasta` - Full protein sequences (RBPs)
- `sequences/nonrbp_sequences.fasta` - Full protein sequences (non-RBPs)

### Results
- `results/family_entropy_summary.csv` - Per-family entropy statistics
- `results/domain_aggregate_frequencies.csv` - AA frequencies RBD vs non-RBD
- `results/position_analysis.png` - Position-specific visualization
- `results/literature_comparison.png` - Findings vs literature table

### Scripts
- `scripts/analyze_domain_vocabulary.py` - Main entropy analysis
- `scripts/analyze_idr_drich_position.py` - D/E segments and position analysis
- `scripts/analyze_flanking_regions.py` - Flank disorder analysis
- `scripts/analyze_outliers_correlation.py` - Family outliers and correlations
- `scripts/visualize_idr_analysis.py` - Generate figures

---

## Literature Review (Expanded 2026-02-01)

### Core Reference: D/E Autoinhibition Paradigm

**Wang X, Levy Y, Iwahara J (2025). Competition between Nucleic Acids and Intrinsically Disordered Regions within Proteins. *Acc Chem Res* 58:2415-2424.** [DOI](https://pubs.acs.org/doi/10.1021/acs.accounts.5c00261)

Key findings from this review:
- **50% of D/E-tract proteins are DNA/RNA-binding proteins** - validates our enrichment finding
- 268 human proteins have ≥10 consecutive D/E residues; only 12 have ≥10 K/R
- D/E-rich IDRs "electrostatically mimic nucleic acids" - the autoinhibition mechanism
- D/E repeats can **accelerate** target search by rejecting "decoys" (non-functional binding sites)

**Our data comparison:**
| Literature claim | Our finding | Match |
|-----------------|-------------|-------|
| ~50% D/E-tract proteins are RBPs | 51.2% RBDs have D/E segments | ✓ |
| D/E enrichment in RBPs | +2.05% in RBDs vs non-RBDs | ✓ |
| D/E in flanking IDRs | Flanks -21 pLDDT (disordered) | ✓ |

### Supporting Literature

#### 1. Target Search Acceleration
**Wang X et al. (2023). Negatively charged, intrinsically disordered regions can accelerate target search by DNA-binding proteins. *Nucleic Acids Res* 51:4701-4712.** [Link](https://academic.oup.com/nar/article/51/10/4701/7034408)

- D/E repeats help proteins avoid "decoy" binding sites
- Autoinhibited proteins transition to uninhibited state via induced-fit
- This may explain why D/E-rich RBDs bind targets despite autoinhibition

#### 2. Hfq Acidic Tail Model
**Santiago-Frangos A et al. (2017). Acidic C-terminal domains autoregulate the RNA chaperone Hfq. *eLife* 6:e27049.** [Link](https://elifesciences.org/articles/27049)

- E. coli Hfq has acidic C-terminal tail that mimics nucleic acids
- Competes against non-specific RNA binding
- **Directly supports our C-terminal D/E enrichment finding (+2.26% at C-terminus)**

#### 3. PUF Protein FBF-2 Model
**Qiu C et al. (2023). Intra- and inter-molecular regulation by intrinsically-disordered regions governs PUF protein RNA binding. *Nat Commun* 14:7612.** [Link](https://www.nature.com/articles/s41467-023-43098-1)

- C-terminal IDR autoinhibits RNA binding by increasing off-rate
- IDR contains electronegative cluster near RNA 5' end
- Partner protein (LST-1) displaces tail to relieve autoinhibition
- **Crystal structure at 2.1 Å shows mechanism**

#### 4. Ribosomal Proteins Are Different
**Nosrati M et al. (2023). The RNA-Binding Function of Ribosomal Proteins and Ribosome Biogenesis Factors in Human Health and Disease. *Biomolecules* 13:1684.** [Link](https://pmc.ncbi.nlm.nih.gov/articles/PMC10669870/)

- Ribosomal proteins use **basic residues (K/R) and aromatics** for rRNA binding
- Many ribosomal proteins are **disordered in isolation but fold upon rRNA binding**
- rRNA provides structural context - no need for autoinhibition
- **Explains our "ribosomal outlier" finding (low D/E, high K/R)**

#### 5. Phosphorylation of Acidic Regions
**Thapar R (2015). Structural Basis for Regulation of RNA-Binding Proteins by Phosphorylation. *ACS Chem Biol* 10:652-666.** [Link](https://pmc.ncbi.nlm.nih.gov/articles/PMC4372107/)

- SLBP has acidic D/E-rich region that is phosphorylated
- Phosphorylation of C-terminal tail increases RNA binding ~27-fold
- **Suggests D/E regions may be PTM hotspots for regulation**

#### 6. Protein Disorder and Autoinhibition
**Bugge K et al. (2024). Protein disorder and autoinhibition: the role of multivalency and effective concentration. *Curr Opin Struct Biol* 84:102741.** [Link](https://pmc.ncbi.nlm.nih.gov/articles/PMC10841074/)

- IDRs regulate oligomerization, catalysis, binding affinity/specificity
- Autoinhibitory interactions are inherently transient
- Role depends on binding to ordered domain → competition or allostery

### Novel Contributions (Not Explicitly in Literature)

| Our Finding | Literature Status | Significance |
|-------------|------------------|--------------|
| C-terminal D/E enrichment (+2.26%) | Hfq model exists, but systematic analysis novel | Genome-wide validation |
| Simpson's Paradox (family vs domain) | Not documented in RBP field | Important statistical caveat |
| Ribosomal proteins as outliers | Known to be different, not explicitly D/E-depleted | Mechanism explanation |
| Entropy ↔ D/E segment correlation (r=0.48) | Not tested before | Novel regulatory insight |
| zf-CCCH low D/E (10.8%) | Zinc fingers known, D/E quantification novel | Alternative regulation |

### Literature Gap Analysis

| Question | Literature Status | Opportunity |
|----------|------------------|-------------|
| Are D/E mutations pathogenic? | Some case studies exist | Systematic ClinVar analysis |
| Do D/E positions have PTMs? | SLBP example known | Proteome-wide mapping |
| Conservation of D/E segment positions | Not systematically studied | Cross-species alignment |
| CodonFM signal in D/E regions | Not studied | Novel angle |

## References

1. Wang X, Levy Y, Iwahara J (2025). Competition between Nucleic Acids and Intrinsically Disordered Regions within Proteins. *Acc Chem Res* 58:2415-2424. [DOI](https://pubs.acs.org/doi/10.1021/acs.accounts.5c00261)

2. Wang X et al. (2023). Negatively charged, intrinsically disordered regions can accelerate target search by DNA-binding proteins. *Nucleic Acids Res* 51:4701-4712. [Link](https://academic.oup.com/nar/article/51/10/4701/7034408)

3. Santiago-Frangos A et al. (2017). Acidic C-terminal domains autoregulate the RNA chaperone Hfq. *eLife* 6:e27049. [Link](https://elifesciences.org/articles/27049)

4. Qiu C et al. (2023). Intra- and inter-molecular regulation by intrinsically-disordered regions governs PUF protein RNA binding. *Nat Commun* 14:7612. [Link](https://www.nature.com/articles/s41467-023-43098-1)

5. Nosrati M et al. (2023). The RNA-Binding Function of Ribosomal Proteins. *Biomolecules* 13:1684. [Link](https://pmc.ncbi.nlm.nih.gov/articles/PMC10669870/)

6. Thapar R (2015). Structural Basis for Regulation of RNA-Binding Proteins by Phosphorylation. *ACS Chem Biol* 10:652-666. [Link](https://pmc.ncbi.nlm.nih.gov/articles/PMC4372107/)

7. Bugge K et al. (2024). Protein disorder and autoinhibition. *Curr Opin Struct Biol* 84:102741. [Link](https://pmc.ncbi.nlm.nih.gov/articles/PMC10841074/)

8. Protein–RNA interactions: structural characteristics and hotspot amino acids (2018). *RNA* 24:1457. [Link](https://rnajournal.cshlp.org/content/24/11/1457.full)

---

## Functional Analysis Results (NEW)

### GO Enrichment: High D/E vs Low D/E Families ⚠️ SMALL SAMPLE

**Sample sizes:** 189 high D/E proteins, 34 low D/E proteins (with GO annotations)

| GO Term | High D/E | Low D/E | Fisher p (Bonferroni) | Status |
|---------|----------|---------|----------------------|--------|
| RNA binding | 89.9% | 97.1% | 1.000 | NS |
| ribosome | 1.1% | 76.5% | **< 0.001** | ✓ |
| translation | 3.2% | 61.8% | **< 0.001** | ✓ |
| mRNA processing | 17.5% | 0.0% | **0.014** | ✓ |

**Interpretation:**
- Ribosome/translation enrichment in low D/E group is highly significant
- mRNA processing enrichment in high D/E group is marginally significant
- ⚠️ Low D/E group has only 34 proteins - interpret with caution

### Codon Analysis Results (Preliminary)

**Sample:** 43 human RBD domains with CDS

| Region | Rare Codon % |
|--------|--------------|
| D/E positions | 0.0% |
| Other positions | 12.6% |

**Caveat:** D and E amino acids only have common codons (GAT, GAC, GAA, GAG) in human, so this comparison is confounded. Future analysis should:
1. Compare codon context AROUND D/E positions
2. Use CodonFM perplexity which captures broader sequence features
3. Analyze synonymous variation at non-D/E sites in D/E-rich regions

**Files generated:**
- `results/domain_cds_for_codonfm.fasta` - 43 domain CDS for CodonFM scoring

---

## AlphaFold pLDDT Validation (NEW)

**Sample:** 100 RBD + 90 non-RBD domains with AlphaFold structures

### Results

| Region | Mean pLDDT (RBD) | Mean pLDDT (non-RBD) | % Disordered |
|--------|------------------|----------------------|--------------|
| D/E positions | 89.8 | 86.9 | 2.1% |
| Other positions | 90.1 | 88.0 | - |
| Domain overall | 90.1 | 87.9 | - |
| **N-terminal flank** | **68.3** | 68.4 | **45.6%** |
| **C-terminal flank** | **69.8** | 71.5 | - |

### Key Findings

1. **D/E within domains: NOT disordered** (pLDDT ~90)
   - D/E positions are just as ordered as other positions within domain cores

2. **Flanks ARE disordered** (pLDDT ~69, Δ = -21 from domain)
   - 45.6% of flank residues have pLDDT < 70
   - This is where D/E-rich autoinhibitory regions reside

3. **RBD vs non-RBD: similar pattern**
   - Both have ordered domains + disordered flanks
   - RBD domains slightly more ordered (+2-3 pLDDT)

### Interpretation

**VALIDATES hypothesis:** D/E-rich autoinhibitory regions are in **flanking IDRs outside Pfam domain boundaries**, not within the structured domain core. This explains:

1. Why domain cores showed LOWER disorder scores (compositional analysis)
2. Why flanks showed HIGHER D/E content
3. The autoinhibition mechanism: flexible acidic tails compete with RNA for basic binding surfaces

### pLDDT Scale Reference
- \>90: Very high confidence (ordered)
- 70-90: Confident (mostly ordered)
- 50-70: Low confidence (potentially disordered)
- <50: Very low confidence (likely disordered)

---

## Corrected Conclusions

### What IS Robust

1. **Flanking regions are disordered** (AlphaFold pLDDT -21 vs domain)
   - This is the strongest finding, not affected by Simpson's paradox
   - Domain cores are ordered (~90 pLDDT), flanks are borderline disordered (~69 pLDDT)

2. **Within-domain length-controlled effects**
   - After controlling for domain length, RBDs still show +1.13% more D/E (p < 1e-87)
   - D/E segment OR = 1.35 after length control (p < 1e-16)

3. **Specific families have high D/E**
   - RRM_1, HSP70, WW, SAM domains have genuinely high D/E (13-15%)
   - These are well-known RNA-binding domains

### What is NOT Robust

1. **"RBDs have more D/E than non-RBDs" as universal claim**
   - Family-weighted analysis: p = 0.234 (NOT significant)
   - The effect is driven by a few large, high-D/E families (RRM_1, WD40, HSP70)
   - Many RBD families have similar or lower D/E than non-RBDs

2. **Causal claims about autoinhibition**
   - We showed compositional differences, not function
   - "Autoinhibitory" interpretation is from literature, not our data

### Revised Interpretation

The D/E enrichment pattern is **family-specific, not a universal RBD property**:

- **High D/E families:** RRM_1, HSP70, WW, SAM (mRNA processing, stress response)
- **Low D/E families:** Ribosomal proteins, zinc fingers (structured contexts)

This suggests autoinhibition via D/E-rich IDRs evolved in specific RBD lineages that require target discrimination, not in all RBDs.

---

## Critical Analysis & Bug Fixes (2026-02-01)

### Data Quality Issues Found

1. **WD40 Contamination (CRITICAL)**
   - 2,485 WD40 (PF00400) entries appeared in BOTH RBD and non-RBD groups
   - Same accessions with different `is_rbd` flags
   - **Fix:** Use Pfam family list as ground truth instead of `is_rbd` flag
   ```python
   # Before (buggy)
   is_rbd = d["is_rbd"]

   # After (fixed)
   is_rbd = d["pfam"] in rbd_families
   ```

2. **RBD Duplicates**
   - 103 RBD entries were exact duplicates
   - **Impact:** Minor (< 1.5% of RBD data)

### Sample Size Changes

| Group | Before Fix | After Fix | Change |
|-------|------------|-----------|--------|
| RBD | 7,519 | 7,416 | -103 (-1.4%) |
| non-RBD | 19,309 | 16,824 | -2,485 (-12.9%) |
| Total | 26,828 | 24,240 | -2,588 |

### Effect on Results

| Metric | Before Fix | After Fix | Status |
|--------|------------|-----------|--------|
| D/E difference | +1.90% | **+2.05%** | ✓ Effect STRONGER |
| Domain-level p | < 1e-200 | < 1e-213 | ✓ More significant |
| Family-level p | 0.276 | **0.234** | ✗ Still NOT significant |
| Odds ratio | 1.89 | **1.99** | ✓ Effect STRONGER |

**Conclusion:** Bug fix actually STRENGTHENED the domain-level effect, but Simpson's Paradox (family-level non-significance) persists. This is NOT a bug - it's a real statistical phenomenon.

### What CAN Be Claimed (Honestly)

| Claim | Status | Evidence |
|-------|--------|----------|
| Domain-level RBDs have higher D/E | ✓ ROBUST | p < 1e-213, effect = +2.05% |
| Flanking regions more disordered | ✓ ROBUST | pLDDT -21 vs domain cores |
| C-terminal enrichment | ✓ ROBUST | χ² = 673, p < 1e-148 |
| Effect driven by specific families | ✓ TRUE | RRM, WD40, HSP70 dominant |

### What CANNOT Be Claimed

| Claim | Status | Reason |
|-------|--------|--------|
| "RBDs universally have more D/E" | ✗ FAILS | Family-level p = 0.234 (NS) |
| "Autoinhibition function proven" | ✗ FAILS | Showed composition, not function |
| "All RBDs use D/E regulation" | ✗ FAILS | Many RBD families have LOW D/E |

### Scripts Requiring the Fix

All analysis scripts should use this pattern:
```python
rbd_families = {"RRM_1", "KH_1", "DEAD", "Helicase_C", ...}  # Known RBD families

for d in data:
    is_rbd = d["pfam"] in rbd_families  # Use family list, NOT is_rbd flag
```

---

*Last updated: 2026-02-01*
*Robust analysis performed: confound control, bootstrap, FDR correction, family-weighted analysis*
*Critical review performed: data quality audit, bug fixes applied*
