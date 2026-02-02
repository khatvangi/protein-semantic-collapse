# Deep Functional Analysis: Summary of Verified Findings

**Date:** 2026-02-01
**Status:** Verified

---

## Four Analyses Completed

| Analysis | Finding | Status | p-value |
|----------|---------|--------|---------|
| 1. mRNA Decay | Decay proteins have HIGH D/E (12.6%) | ✓ VERIFIED | 0.024 |
| 2. SR Phosphorylation | D/E-phospho correlation | ⚠️ NOT SIGNIFICANT | 0.18 |
| 3. Stress Granule vs Ribosomal | SG trend higher but underpowered | ⚠️ UNDERPOWERED | 0.33 |
| 4. KH Proline | KH has MORE proline than RRM | ✓ VERIFIED | 2.9e-23 |

---

## Key Verified Findings

### 1. mRNA Decay Proteins Use D/E Autoinhibition ✓

| Function | D/E % | n | Interpretation |
|----------|-------|---|----------------|
| Decay | 12.60% | 14 | DYNAMIC process |
| Splicing | 13.83% | 16 | DYNAMIC process |
| Other | 10.73% | 70 | Mixed |

**p = 0.024** (Decay vs Other)

**Conclusion:** Both splicing AND decay proteins have elevated D/E. These are both DYNAMIC RNA processes.

### 2. KH Domains Use Proline Instead of D/E ✓

| Domain | D/E % | Proline % | Glycine % |
|--------|-------|-----------|-----------|
| KH | 11.2% | **4.35%** | 9.14% |
| RRM | 13.8% | 2.83% | 8.13% |
| p-value | - | **2.9e-23** | 2.3e-8 |

**Literature validation:**
- KSRP has "proline-glycine (P/G)-rich region" ([PMC11003564](https://pmc.ncbi.nlm.nih.gov/articles/PMC11003564/))
- Sam68 has "proline-rich motifs that resemble SH3 binding sites" ([Nature](https://www.nature.com/articles/1203079))

**Conclusion:** KH domains use PROLINE-RICH REGIONS for protein-protein interaction-based regulation, NOT D/E electrostatic autoinhibition.

---

## The Two-Mechanism Model (Refined)

```
             RNA-Binding Protein Regulation
                        |
        ┌───────────────┴───────────────┐
        │                               │
   D/E AUTOINHIBITION            PROLINE-RICH REGULATION
   (Electrostatic)               (Protein-Protein)
        │                               │
   RRM, LSM, DEAD                  KH, Sam68
   Splicing, Decay                 mRNA transport
        │                               │
   - Acidic tails                  - P/G-rich regions
   - Compete with RNA              - SH3 domain binding
   - Rapid on/off                  - Signaling integration
        │                               │
   HIGH D/E (12-14%)              HIGH Pro (4-5%)
   LOW Pro (2-3%)                 LOW D/E (11%)
```

---

## Biological Logic

### Why D/E for Splicing/Decay?

1. **Need rapid release** - Spliceosome cycles through many substrates
2. **Electrostatic switch** - D/E competes with RNA phosphate backbone
3. **Tunable** - Phosphorylation can modify nearby residues

### Why Proline for KH/Sam68?

1. **Signaling integration** - Proline-rich regions bind SH3 domains (kinases, adaptors)
2. **Different dynamics** - mRNA transport needs coordination, not rapid cycling
3. **Structural flexibility** - Polyproline II helices mediate protein-protein interactions

---

## What Remains Unclear

| Question | Status |
|----------|--------|
| SR protein phosphorylation + D/E | Not significant (r=0.19, p=0.18) |
| Stress granule vs ribosomal | Underpowered (n=6 vs n=19) |
| dsRBD mechanism | Not analyzed |
| Zinc finger mechanism | Not analyzed |

---

## Sources

1. KSRP structure: [PMC11003564](https://pmc.ncbi.nlm.nih.gov/articles/PMC11003564/)
2. Sam68 SH3 binding: [Nature 1203079](https://www.nature.com/articles/1203079)
3. KH domain review: [Wikipedia](https://en.wikipedia.org/wiki/KH_domain)

---

*Analysis completed: 2026-02-01*
