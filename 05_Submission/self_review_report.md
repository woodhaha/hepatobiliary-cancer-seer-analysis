# Self-Review Report: Comparable but Not Equivalent — Segmental vs Major Hepatectomy in Elderly Hepatobiliary Cancer

**Target journal**: Annals of Surgical Oncology
**Manuscript type**: Observational cohort (STROBE + TRIPOD)
**Date**: 2026-06-26
**Overall assessment**: Methodologically robust (E-value=7.5, 76K patients, dual external validation). Three fixable Major issues before submission: (1) DeepSurv model failure unexplained, (2) ICGC external C-index near random, (3) unsubstantiated informed-consent statement. Minors are mostly reporting-compliance items.

---

## Deterministic Gate Results

| Gate | Result |
|------|--------|
| Scope Coherence (§D) | ✅ PASS — conclusion matches design/endpoint |
| Classical Style (§J) | ✅ PASS — no AI tells, § marks, or policy violations |
| Cohort Arithmetic (Phase 2.5) | ✅ PASS — STROBE flow balances, no n-mismatch |

---

## Anticipated Major Comments (fix before submission)

**M1. DeepSurv model failure is unexplained and undermines ML claims** [C. Validation]
The manuscript reports DeepSurv C-index = 0.575, barely above random (0.5), while all other models achieve 0.67–0.75. This is a 12–18 point performance drop — far too large for normal model variance. A reviewer will suspect: (a) implementation error in the Cox partial-likelihood loss, (b) insufficient tuning (the model achieved only 0.523 before tuning), or (c) data leakage in the tree-based models inflating their performance relative to the neural network. The manuscript currently reports DeepSurv without comment. Either fix the implementation, document extensive tuning attempts with negative results, or remove DeepSurv and reframe as "3-model ensemble" rather than "5-model comparison."

**Severity**: Fixable
**Suggested fix**: (a) Remove DeepSurv from the main comparison and move to Supplementary as an exploratory analysis with explicit limitations noted. OR (b) Re-implement with proper learning-rate schedule, wider architecture search, and document convergence diagnostics. The current 0.575 cannot stand as-is in a JAMA Surgery submission.
**fixable_by_ai**: false (requires re-running analysis)

**M2. ICGC external C-index = 0.522 (Cox) is at chance level — need honest framing** [C. Validation]
The ICGC-LIRI-JP validation yields C-index 0.522 for the Cox model and 0.547 for XGBoost. These are near or at chance level (95% CI likely crosses 0.5). The manuscript acknowledges the gap but attributes it entirely to the "Surgery Feature Paradox." A reviewer will note that TCGA (also all-surgical) achieved 0.595 — a 7-point gap between the two external cohorts that the surgery-paradox explanation cannot fully explain. The ICGC cohort's Japanese population introduces population-drift (different etiology, baseline risk, staging practices) that the SEER model does not capture. This should be explicitly discussed as a population-generalizability limitation, not solely a feature-variance issue.

**Severity**: Fixable
**Suggested fix**: Add a paragraph in Discussion explicitly addressing the TCGA-ICGC differential: "The 7-point C-index gap between TCGA (0.595) and ICGC (0.522) cannot be attributed solely to the Surgery Feature Paradox, as both are all-surgical. Population differences—Japanese vs US etiology (HBV-dominant vs HCV/alcohol), staging calibration, and baseline risk profiles—likely contribute. ICGC performance near chance underscores that SEER-trained models require population-specific recalibration before Asian-population deployment."
**fixable_by_ai**: true

**M3. "Informed consent was obtained" statement is unsubstantiated** [F. Reporting]
The abstract states "Informed consent was obtained from all patients in the external clinical validation cohort" — but this study has no Wenzhou clinical validation cohort. This appears to be a residual from the gastric cancer project template. A reviewer or editor will flag this immediately as a factual error in an ethics statement. If a clinical cohort exists and was used, it must be described in Methods with N, dates, and IRB approval. If not, the statement must be removed.

**Severity**: Fixable (but P0 before submission — ethics statement error)
**Suggested fix**: Remove the informed-consent sentence. SEER, TCGA, and ICGC are all publicly available de-identified databases that do not require individual patient consent. Replace with: "SEER, TCGA, and ICGC data are publicly available and de-identified; institutional review board approval was exempt."
**fixable_by_ai**: true

**M4. Missing TRIPOD compliance for prediction model development** [G. Reporting Guideline]
The study develops a clinical prediction model (HBI risk score, nomogram, ensemble survival model). TRIPOD (Transparent Reporting of a multivariable prediction model for Individual Prognosis Or Diagnosis) compliance is mandatory for journals at this tier. Key TRIPOD items not addressed: (a) blinding, (b) model specification rationale, (c) handling of predictors, (d) model presentation for individual predictions. While complete TRIPOD+AI may be excessive (the ML component is secondary), basic TRIPOD items should be addressed.

**Severity**: Fixable
**Suggested fix**: Add a TRIPOD compliance statement in Methods: "This study follows the TRIPOD (Transparent Reporting of a multivariable prediction model for Individual Prognosis Or Diagnosis) statement for prediction model development and validation (Type 3: development + external validation)." Address items 13a (model presentation), 15a (full model specification), and 16 (model usage).
**fixable_by_ai**: true

**M5. IV-adjusted HR = 0.05 is implausible and suggests weak instrument** [A. Study Design]
The instrumental variable analysis using regional income as instrument yields an implausibly extreme HR of 0.05 (vs PSM-adjusted 0.35 and multivariate 0.25). This 5–7× stronger effect violates the monotonicity assumption and indicates a weak instrument. Reporting this without qualification invites dismissal of the entire IV analysis. The income instrument is inherently weak: while income correlates with surgery receipt (R²=0.06), it also directly affects survival through mechanisms independent of surgery (access to care, nutrition, follow-up compliance).

**Severity**: Fixable
**Suggested fix**: Either (a) remove the IV analysis and replace with a more robust causal method (e.g., inverse probability of treatment weighting with stabilized weights), or (b) explicitly flag the IV as exploratory: "Given the weak instrument (Stage 1 pseudo-R²=0.06), the IV estimate should not be interpreted as the true causal effect but as directional confirmation of the PSM and E-value findings." The current unqualified reporting of HR=0.05 is misleading.
**fixable_by_ai**: true (text-only); requires_reanalysis if switching to IPTW

**M6. No calibration assessment for prediction models** [C. Validation]
The manuscript reports discrimination (C-index) for all models but no calibration. For prediction models intended for clinical use (nomogram, HBI risk score), calibration is equally or more important than discrimination. JAMA Surgery reviewers routinely require calibration plots, calibration-in-the-large, and calibration slope. The bootstrap calibration framework is mentioned in the analysis outputs but not reported in the manuscript.

**Severity**: Fixable
**Suggested fix**: Add calibration section to Results reporting: (1) calibration-in-the-large, (2) calibration slope, (3) Brier score at 12/36/60 months for the ensemble model. Include a calibration plot (predicted vs observed survival) as a supplementary figure. If calibration is poor despite good discrimination, this is important to discuss honestly.
**fixable_by_ai**: false (requires generating calibration figures from data)

---

## Anticipated Minor Comments (address proactively)

**m1. "Other" category (6,981, 9.2%) never defined** [F. Reporting]: The cohort includes 6,981 "other" patients neither HCC nor ICC, but they are never described. Specify: what are these patients? (e.g., combined HCC-ICC, hepatoblastoma, sarcoma NOS). Either define them or, if they were excluded from the primary analysis, update the CONSORT diagram accordingly.

**m2. Gradient Boosting Survival mentioned but not run** [E. Reproducibility]: Methods lists 5 models including Gradient Boosting, but the analysis actually skipped GB due to computational constraints (~870s training time). Either run it (even on a subsample) or remove from the model list and state "4 models were trained."

**m3. Missing STROBE checklist** [G. Reporting]: Observational study should include a completed STROBE checklist as supplementary material. This is standard for surgical journals.

**m4. Median follow-up with IQR not reported** [F. Reporting]: The manuscript reports median follow-up (49m OS, 74m CSS) but without IQR. For survival studies, median follow-up should be reported with IQR using the reverse Kaplan-Meier method.

**m5. The "paradigm shift" language is strong for observational data** [D. Clinical Framing]: The phrase "paradigm shift" appears in Conclusions. While the evidence is robust (E-value=7.5), "paradigm shift" is editorializing. JAMA Surgery reviewers prefer more measured language. Suggest: "These findings support parenchymal-sparing surgery as the preferred approach."

**m6. SEER covers 18 registries, not all reporting to ~28%** [A. Study Design]: Methods says "SEER covers approximately 28% of the US population across 18 registries" — this is accurate for the 2024 submission. The introduction paraphrases this correctly. Verify the exact coverage percentage for the November 2024 submission.

**m7. No data on surgical approach (open vs laparoscopic)** [D. Clinical]: SEER does not distinguish laparoscopic from open resection. This is a significant limitation for a surgical journal — laparoscopic segmentectomy outcomes may differ from open. Add to Limitations.

**m8. No assessment of lymph node dissection extent** [D. Clinical]: For ICC, adequate lymphadenectomy (≥6 nodes) is a known quality metric. SEER has nodes examined/positive data. Add a brief note on whether this was available and why it was not modeled.

---

## Strengths (emphasize in cover letter)

- **Exceptional sample size**: 76,110 elderly patients — 32× larger than Zhang et al. (2020), 7.5× larger than Xia et al. (2025)
- **Causal rigor unprecedented in surgical SEER studies**: E-value=7.5 + PSM (SMD 0.190→0.030) + 6 sensitivity analyses
- **First formal surgery×cancer-type interaction test in hepatobiliary surgery literature**: p=0.206, establishing equivalence across HCC and ICC
- **Dual external validation on real TCGA-LIHC data + ICGC-LIRI-JP**: first SEER hepatobiliary model with external benchmarks
- **Clinical decision matrix**: directly actionable — specifies segmental vs larger by age/cancer/stage/frailty

---

## R0 Pre-Submission Findings (for /revise cross-reference)

R0-1 [MAJ] DeepSurv model failure unexplained — remove or document
R0-2 [MAJ] ICGC near-chance performance needs population-drift discussion
R0-3 [MAJ] Unsubstantiated informed-consent statement — remove
R0-4 [MAJ] Missing TRIPOD compliance statement
R0-5 [MAJ] IV HR=0.05 implausible — weak instrument, flag or remove
R0-6 [MAJ] No calibration assessment for prediction models
R0-7 [MIN] "Other" cancer category undefined
R0-8 [MIN] Gradient Boosting listed but not run
R0-9 [MIN] Missing STROBE checklist
R0-10 [MIN] Median follow-up needs IQR
R0-11 [MIN] "Paradigm shift" too strong — soften language
R0-12 [MIN] Add laparoscopic vs open limitation
R0-13 [MIN] Lymph node dissection for ICC not addressed
