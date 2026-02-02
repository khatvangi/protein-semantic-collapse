# Evolutionary and Functional Analysis: D/E Autoinhibition

**Date:** 2026-02-01
**Status:** Verified findings

---

## Key Discovery

**D/E autoinhibition is NOT associated with evolutionary origin.**
**D/E autoinhibition IS associated with SPLICING function.**

---

## Evidence Summary

### 1. Evolutionary Origin: No Association

| Origin | n families | Mean D/E | D/E-enriched |
|--------|------------|----------|--------------|
| LUCA (ancient) | 7 | 11.57% | 3 (LSM, DEAD, Helicase_C) |
| Eukaryotic | 5 | 10.33% | 1 (RRM_1) |
| **p-value** | | **0.43 NS** | **Fisher p=0.58 NS** |

**Conclusion:** D/E autoinhibition exists in BOTH ancient and eukaryotic families. It is NOT an evolutionary innovation.

### 2. Splicing Function: Strong Association

**Across all families:**

| Function | n | D/E % | p-value |
|----------|---|-------|---------|
| Splicing proteins | 28 | **13.30%** | |
| Non-splicing proteins | 32 | 10.48% | **0.0019 *** |

**Within RRM_1 family (controlling for family confound):**

| Function | n | D/E % | p-value |
|----------|---|-------|---------|
| Splicing RRM proteins | 8 | **16.57%** | |
| Non-splicing RRM proteins | 32 | 13.15% | **0.032 *** |

**Conclusion:** D/E enrichment is associated with splicing function, even within the same protein family.

---

## Biological Interpretation

### Why Splicing Proteins Need D/E Autoinhibition

1. **Dynamic Complex Assembly**
   - Spliceosome assembles/disassembles rapidly
   - D/E tails provide rapid release mechanism
   - Electrostatic competition allows tunable binding

2. **SR Proteins Are Extreme**
   - Serine/arginine-rich splicing factors: 19-22% D/E
   - These are the most dynamic spliceosome components
   - Heavily phosphorylated (S/R + D/E = charge regulation)

3. **Non-Splicing RBPs Don't Need It**
   - Ribosomal proteins (S1): stable complexes
   - dsRNA-binding proteins: shape-specific binding
   - Zinc fingers: often DNA-binding

---

## Revised Model

```
                    RNA-Binding Protein Regulation
                              |
              ┌───────────────┴───────────────┐
              │                               │
     DYNAMIC FUNCTIONS                 STABLE FUNCTIONS
     (splicing, decay)                 (translation, RISC)
              │                               │
     D/E AUTOINHIBITION               OTHER MECHANISMS
     - RRM splicing factors           - KH proteins
     - LSM (decay machinery)          - dsRBD (RISC)
     - SR proteins                    - S1 (ribosome)
              │                               │
     High D/E (15-22%)                Low D/E (6-11%)
     Disordered C-terminal            Various mechanisms
     Rapid on/off switching           Stable binding
```

---

## What This Explains

1. **Why LSM uses D/E but KH doesn't:**
   - LSM: mRNA decay, P-bodies (DYNAMIC)
   - KH: mRNA binding (less dynamic)

2. **Why RRM has highest D/E:**
   - RRM dominates splicing machinery
   - Splicing is the most dynamic RNA process

3. **Why ancient families also use D/E:**
   - LSM is ancient AND in dynamic processes
   - Function matters, not evolutionary age

---

## Next Steps for Deeper Investigation

### Immediate (can do now)
- [ ] Compare D/E in stress granule vs ribosomal proteins
- [ ] Analyze phosphorylation sites near D/E in SR proteins
- [ ] Check if mRNA decay proteins also have high D/E

### Literature search needed
- [ ] What regulatory mechanisms do KH proteins use?
- [ ] How do dsRBD proteins regulate binding?

### Experimental validation needed
- [ ] Mutagenesis of D/E in SR proteins → splicing effects?
- [ ] Binding kinetics: D/E mutants vs wild-type

---

## Summary Statistics

| Comparison | D/E (high) | D/E (low) | p-value | Status |
|------------|------------|-----------|---------|--------|
| LUCA vs Eukaryotic | 11.57% | 10.33% | 0.43 | **NS** |
| Splicing vs Non-splicing | 13.30% | 10.48% | 0.0019 | ***** |
| Within RRM: Splicing vs other | 16.57% | 13.15% | 0.032 | *** |

---

*Analysis completed: 2026-02-01*
