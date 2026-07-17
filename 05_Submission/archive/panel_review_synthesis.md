# Panel Review Synthesis — Hepatobiliary Cancer Manuscript

**Panel composition**: R1 (Epidemiology/Statistics) · R2 (Clinical Hepatobiliary Surgery) · R3 (ML/AI Methods)
**Date**: 2026-06-26
**Editorial Verdict**: **REVISE** — All 3 reviewers identified fixable issues. No fatal design flaws.

---

## CRITICAL / Must-Fix (P0) — 3 items

### C1. "Corrected Cox partial likelihood" is a fabricated term [R3 · CONSENSUS R1+R3]
Both R1 and R3 identified the DeepSurv loss description as problematic. The standard term is "negative log partial likelihood of the Cox proportional hazards model." The word "corrected" implies a modification that does not exist in the literature. R1's code inspection confirms the bug was fixed (sort by time, not risk) but the term is still wrong.
**Fix**: Replace with "negative log partial likelihood of the Cox model." Document the sort-by-time fix in supplementary methods.

### C2. Reference 27 is completely unrelated [R2, R3]
Ormeño et al. (Commun Chem, 2024) is about molecular dynamics simulations. Zero relevance to surgical oncology.
**Fix**: Remove ref 27. Verify ALL 28 references against DOIs.

### C3. E-value=7.5 is overstated — computed on wrong HR [R1, R2]
The E-value was computed on the multivariate HR of 0.23 (surgery_any in a full model including mediators). The appropriate conservative estimate is on the PSM-adjusted HR of 0.35, yielding **E-value ≈ 4.3**. R2 further notes the E-value addresses surgery-vs-none, not the segmental-vs-larger comparison that is the paper's main claim.
**Fix**: Report E-value=5.2 (PSM) as primary. Keep E-value=7.5 (multivariate) as secondary. Add E-value for the segmental-vs-larger specific comparison.

---

## MAJOR — 5 items

### M1. Fine-Gray analysis is actually cause-specific hazard [R1]
The code uses `CoxPHFitter` on filtered subsets — cause-specific hazards, not Fine-Gray subdistribution.
**Fix**: Rename "Fine-Gray" → "Cause-specific hazard regression." Implement true Fine-Gray OR explain why CSH is preferred.

### M2. Confounding by indication not addressed for the segmental-vs-larger comparison [R2]
The E-value covers surgery-vs-none. The segmental-vs-larger equivalence may reflect appropriate surgical selection, not oncologic equivalence.
**Fix**: Add Discussion paragraph acknowledging this. Soften "default strategy" → "when both approaches are anatomically feasible."

### M3. Ensemble excludes DeepSurv (best model) without justification [R3]
DeepSurv C=0.750 but the ensemble combines only Cox+RSF+XGBoost, yielding C=0.737 < all constituents.
**Fix**: Include DeepSurv in ensemble OR justify exclusion. If C stays <0.750, discuss why.

### M4. Liver function variables absent [R2]
SEER lacks Child-Pugh/MELD — the dominant confounder in HCC surgery.
**Fix**: Add dedicated paragraph to Discussion. The cirrhosis paradox needs its own explanation.

### M5. External validation ICGC = chance performance [R1, R2]
C-index 0.522-0.547 is at chance level. The manuscript frames this too generously.
**Fix**: Add 95% CIs for external C-indices. Explicitly state "models failed to validate in ICGC." Clarify if ICGC data is synthetic or real individual-level.

---

## MINOR — 8 items

| # | Reviewer | Issue | Fix |
|---|----------|-------|-----|
| m1 | R1 | PH violations unreported | Add Schoenfeld table to supplementary |
| m2 | R3 | Hyperparameter reporting insufficient | Add full table (batch_size, epochs, widths, seeds) |
| m3 | R2 | Surgical categorization collapses wedge+segmental | Discuss limitation; if codes permit, separate |
| m4 | R2 | HBI not clinically actionable as presented | Provide formula + worked example |
| m5 | R3 | "Surgery Feature Paradox" — not a paradox | Rename to "predictor range restriction" |
| m6 | R1 | Immortal time bias unacknowledged | Add to Limitations |
| m7 | R3 | Abstract cohort math: 57,380+11,749≠76,110 | Add "6,981 other" |
| m8 | R2 | BCLC guidelines engagement superficial | Add BCLC-contextualized paragraph |

---

## STRENGTHS (Consensus)

- Impressive sample size (N=76,110) — R1, R2, R3 all note this as exceptional
- E-value framework is a genuine methodological strength — R1, R2
- Clinical question is important and timely — R2
- 5-fold CV + temporal validation design is sound — R3
- ICC inclusion is a genuine novelty — R2
- Calibration assessment is thorough — R1

---

## READINESS VERDICT

**Not ready for Annals of Surgical Oncology in current form.** Fix C1-C3 before re-circulation. Address M1-M5 for full submission readiness. The manuscript's core value (large SEER analysis of elderly hepatobiliary surgical outcomes) is strong; the fixes are all text-level or re-analysis of existing data.
