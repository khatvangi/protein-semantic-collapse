# D/E Autoinhibition Project: Final Status

**Date:** 2026-02-01
**Status:** COMPLETE

---

## Final Results

| Task | Status | Key Finding |
|------|--------|-------------|
| **CodonFM** | ✓ SUCCESS | D/E codons MORE optimized (p≈0) |
| **Write-up** | ✓ SUCCESS | Formal report generated |
| **Structural** | ✓ SUCCESS | C-terminal D/E disordered (pLDDT=36.6) |

---

## Summary of Evidence

### Strong Support ✓

| Finding | Evidence |
|---------|----------|
| D/E removal pathogenic | 1.35x (p=0.007) |
| Amidation pathogenic | 2.3x (p=3.2×10⁻⁵) |
| C-terminal enrichment | +2.26% (p<10⁻¹⁴⁸) |
| C-terminal D/E disordered | pLDDT=36.6 (hnRNPA1) |
| D/E-K/R contacts | 38 within 10Å |
| D/E codons optimized | p≈0 (n=1021) |

### Critical Revision ⚠️

| Original | Revised |
|----------|---------|
| D/E segments MORE common in RBDs (51%) | LESS common (6.7%, OR=0.83) |
| Universal RBD property | Family-specific |

---

## Key Conclusions

1. **Autoinhibition uses DISPERSED D/E**, not consecutive runs
2. **C-terminal D/E regions are disordered** (structural validation)
3. **D/E removal is pathogenic** (functional validation)
4. **D/E codons are optimized** (not under relaxed selection)
5. **Family-specific**: RRM, LSM, DEAD families, not universal

---

## Files Generated

```
results/
├── formal_report.md              # Publication-ready summary
├── strengthened_analysis.txt     # Fixed D/E detection
├── clinvar_de_analysis.txt       # Pathogenicity analysis
├── ptm_de_analysis.txt           # PTM enrichment
├── conservation_de_analysis.txt  # Cross-species conservation
├── mutation_predictions.txt      # D/E→A effect predictions
├── structural_analysis_hnRNPA1.txt  # AlphaFold analysis
└── final_status.json             # Machine-readable summary

structures/
├── hnRNPA1_AF.pdb               # AlphaFold structure
└── PTBP1_AF.pdb                 # AlphaFold structure

HYPOTHESIS_REVIEW.md              # Critical hypothesis revision
FINDINGS.md                       # Complete findings documentation
```

---

## What Would Further Strengthen

1. **MD simulation** (2-5 days GPU time) - Visualize D/E-RNA competition
2. **Wet lab mutagenesis** - Direct binding assays with D/E→A mutants
3. **CLIP-seq integration** - Map D/E to actual RNA binding sites

---

*Project completed: 2026-02-01*
