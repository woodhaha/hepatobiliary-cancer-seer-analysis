# Panel Review Synthesis V2 — Hepatobiliary Cancer Manuscript

**Panel composition**: R1 (Epidemiology/Statistics) · R2 (Clinical Hepatobiliary Surgery) · R3 (ML/AI Methods)
**Date**: 2026-07-11 · **Round**: Post-fix re-review

---

## Editorial Verdict: **ACCEPT WITH MINOR REVISIONS**

The previous round identified 3 critical (P0) and 5 major issues, all of which
have been addressed. The manuscript is substantially improved and suitable for
publication in Annals of Surgical Oncology pending the minor items below.

---

## R1 — Epidemiology & Statistics

### Previous critical issues: RESOLVED

| Issue | Verdict | Evidence |
|-------|---------|----------|
| C3: E-value overstated (7.5→5.2) | ✅ Fixed | E-value 5.2 (PSM) primary, 7.5 (MV) secondary |
| M1: Fine-Gray misnamed | ✅ Fixed | Correctly labeled "cause-specific hazard" |
| M6: No calibration | ✅ Fixed | Brier scores 12/36/60m reported |
| m1: PH violations unreported | ✅ Fixed | Schoenfeld summary in Results + eTable 7 |

### New observations

**1. eFigure 9 referenced but supplementary lists only eFigure 1-8** [Minor]
The manuscript body cites "eFigure 9" for model calibration (Brier scores), but the
Supplementary Materials section lists eFigure 1-8 only. Fig25_Calibration exists
(6 panels: risk stratification, calibration plot, C-index, Brier scores) — either add
it as eFigure 9 or renumber to eFigure 8 and update the "1-8" range.

**2. Power analysis is a welcome addition** [Positive]
The added power paragraph (Methods §Statistical Power) correctly notes >99% power
for the primary comparison and 80% for the surgical subgroup comparison. The caveat
about lower power in individual age-cancer strata is appropriate.

**3. Missing-data handling could be more explicit** [Minor]
The Methods state "complete data on key covariates" as an inclusion criterion, which
implies complete-case analysis. Given the large sample (N=76,110) this is defensible,
but a brief mention in Limitations would strengthen transparency — especially since
the SEER "cirrhosis" field has >90% missingness (acknowledged in Discussion but not
in Methods).

**4. "Formal interaction test" and "marginal advantage" usage** [Style]
The Results appropriately present the interaction test (p=0.206) and describe
age-stratified HRs. At age 65-69, the difference (0.33 vs 0.35) is described as
"marginal advantage" — appropriate language.

---

## R2 — Clinical Hepatobiliary Surgery

### Previous critical issues: RESOLVED

| Issue | Verdict | Evidence |
|-------|---------|----------|
| C2: Ref 27 unrelated (Ormeño → Marrero) | ✅ Fixed | Now Marrero et al., AASLD 2018 — correct |
| M2: Segmental-vs-larger confounding not addressed | ✅ Fixed | E-value ~1.3 + "cautious interpretation" paragraph |
| M4: Liver function absent | ✅ Fixed | Dedicated paragraph + cirrhosis paradox explained |
| m6: Immortal time bias | ✅ Fixed | Acknowledged in Limitations |
| m8: BCLC engagement shallow | ✅ Fixed | Integrated in Introduction + Discussion |

### New observations

**1. Segmental-vs-larger E-value (1.3) is honestly reported** [Positive]
The addition "a weak unmeasured confounder (RR ~1.3) could nullify the observed
equivalence" is appropriately cautious. This paragraph was the key gap in V1.

**2. Clinical decision matrix is actionable** [Strength]
Table 3 provides explicit recommendations by age/cancer/stage/frailty — this is what
clinicians want to see. The recommendation for ICC (segmental preferred over larger)
is novel and appropriately hedged.

**3. Missing: frailty assessment granularity** [Minor]
The FSI uses only age, stage, and grade — these are available in SEER, but a
true frailty index requires comorbidities (Charlson), functional status (ECOG),
and nutrition (albumin). This is correctly listed as a limitation, but the
paragraph could add that even these 3 variables produce clinically meaningful
stratification (13-month OS benefit in the "frail" group).

**4. Missing: number of lymph nodes harvested as quality metric** [Minor]
For ICC, adequate lymphadenectomy (≥6 nodes) is a quality metric [Ref 23].
The manuscript mentions this in Limitations but doesn't report the median node
count from SEER. If available, a brief note (e.g., "median 3 nodes examined,
IQR 1-6") would strengthen the ICC analysis. If unavailable, say so.

---

## R3 — ML/AI Methods

### Previous critical issues: RESOLVED

| Issue | Verdict | Evidence |
|-------|---------|----------|
| C1: "Corrected Cox partial likelihood" fabricated term | ✅ Fixed | Now "negative log partial likelihood" — standard term |
| M3: Ensemble excluded DeepSurv | ✅ Fixed | 4-model ensemble now includes DeepSurv (C=0.751) |
| M5: ICGC framing too generous | ✅ Fixed | 95% CIs reported; "did not successfully validate" |
| m2: Hyperparameters insufficient | ✅ Fixed | Full parameter list for all 4 models |
| m5: "Surgery Feature Paradox" renamed | ✅ Fixed | "Predictor range restriction" — correct |

### New observations

**1. DeepSurv now competitive (0.751) — fix confirmed** [Positive]
After the sort-by-time correction, DeepSurv C-index rose from 0.575 to 0.751, now
the second-best model. The ensemble achieves 0.756 ± 0.003. No further ML concerns.

**2. Brier scores reported but needs calibration plot reference** [Minor]
The Brier scores (RSF 0.157/0.092/0.059 at 12/36/60mo) are good — note that RSF
calibration at 60mo (0.059) approaches the ideal (Brier=0 for perfect, 0.25 for
useless). However, the calibration plot (Figure 25 Panel D) is cited as "eFigure 9"
which is outside the eFigure 1-8 range declared in Supplementary. Reconcile numbering.

**3. No confidence intervals for external validation C-indices** [Minor]
TCGA C-index 0.595 (95% CI 0.536–0.654) is reported — good. But some external
CIs (RSF, XGBoost) are not individually reported. Standard practice: report all.

**4. DeepSurv architecture: 3-layer MLP is appropriate** [Positive]
17→64→32→16→1 with batch norm and dropout=0.3 is a reasonable architecture for
~76K samples. The training details (Adam, lr=1e-3, patience=15, 200 epochs) are
fully specified — reproducible.

---

## Verdict Summary

### What changed from V1 → V2

| Domain | V1 | V2 |
|--------|------|------|
| Critical (P0) | 3 | **0** |
| Major | 5 | **0** |
| Minor | 8 | **5** |
| Verdict | Not ready | **Accept (minor revisions)** |

### Remaining 5 minor items (pre-submission)

| # | Reviewer | Issue | Fix |
|---|----------|-------|-----|
| m1 | R1 | eFigure 9 outside declared 1-8 range | Renumber or re-list |
| m2 | R1 | Missing-data handling not explicitly in Limitations | Add 1 sentence to Limitations |
| m3 | R2 | Node count for ICC not reported | Add median (IQR) if available |
| m4 | R2 | FSI paragraph could mention 13-month benefit | Add 1 sentence |
| m5 | R3 | External validation CIs: report all models | Add RSF/XGBoost CI to text |

### Strengths re-affirmed

- N=76,110 remains the largest such analysis
- Surgery-vs-none causal robustness (E-value 5.2) is publication-grade
- ICC finding (segmental HR 0.22 vs larger HR 0.26) is genuinely novel
- DeepSurv fix validated → 4-model ensemble credible
- Clinical decision matrix is directly actionable

### Recommended Journal

**Annals of Surgical Oncology** — the manuscript's strengths (surgical outcomes,
population-based, geriatric focus) align well with ASO readership.
JAMA Surgery would require shortening and additional mechanistic depth not
possible with SEER data.
