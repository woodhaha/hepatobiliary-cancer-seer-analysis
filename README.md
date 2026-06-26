# Hepatobiliary Cancer SEER Analysis — Complete Pipeline

**Manuscript**: Surgical Resection and Extent in Elderly Patients With Hepatobiliary Cancer  
**Target Journal**: Annals of Surgical Oncology  
**Data**: SEER Research Plus (Nov 2024) + TCGA-LIHC + ICGC-LIRI-JP  
**Date**: 2026-06-26

---

## Analysis Pipeline (16 scripts)

### Phase 1: Data Preparation
| # | Script | Purpose | Output |
|---|--------|---------|--------|
| 1 | `01_data_exploration.py` | Variable distributions, site/histology/year/age/surgery/stage overview | `01_exploration_report.txt` |
| 2 | `02_data_cleaning.py` | Inclusion/exclusion, feature engineering, unified AJCC staging, surgery classification | `hepatobiliary_elderly_clean.csv` |

### Phase 2: Survival & Cox
| # | Script | Purpose | Output |
|---|--------|---------|--------|
| 3 | `03_survival_analysis.py` | Table 1 baseline, 6-panel KM (OS+CSS by surgery/stage/age/cancer type), univariate+multivariate Cox, stratification by age/stage/cancer type, forest plot | `Fig1_KM`, `Fig2_Forest`, `02_survival_report.md` |

### Phase 3: Machine Learning
| # | Script | Purpose | Output |
|---|--------|---------|--------|
| 4 | `04_ml_pipeline.py` | Cox/RSF/XGBoost/GB/DeepSurv training, temporal split, feature importance | `03_ml_report.md`, `ml_data.npz` |

### Phase 4: Advanced Analyses
| # | Script | Purpose | Output |
|---|--------|---------|--------|
| 5 | `05_advanced_analyses.py` | SHAP, competing risk CIF, RMST, age-surgery spline, temporal trends, cirrhosis paradox, HCC vs ICC, interaction tests | `FigS1-S4`, `04_advanced_report.md` |
| 6 | `06_fetch_external_data.py` | Download TCGA-LIHC (UCSC Xena) + ICGC summary data | `external_validation_cohort.csv` |
| 7 | `07_external_validation.py` | Apply SEER models to TCGA+ICGC, surgery paradox diagnosis, risk stratification | `Fig5_ExternalValidation`, `05_external_validation_report.md` |

### Phase 5: Composite & PSM
| # | Script | Purpose | Output |
|---|--------|---------|--------|
| 8 | `08_final_analysis.py` | 5-fold CV ensemble, PSM (1:1 nearest-neighbor), HBI risk score, nomogram | `Fig6_CompositeAnalysis`, `06_manuscript_summary.md` |

### Phase 6: High-Priority Methods
| # | Script | Purpose | Output |
|---|--------|---------|--------|
| 9 | `09_high_priority.py` | LASSO variable selection, DCA, time-dependent AUC, nomogram, DeepSurv tuning | `Fig7_LASSO`, `Fig8_DCA`, `Fig10_Nomogram` |

### Phase 7: Depth Analyses
| # | Script | Purpose | Output |
|---|--------|---------|--------|
| 10 | `10_depth_analyses.py` | Schoenfeld residuals, subgroup interaction forest, landmark survival, AFP analysis, race disparities, COVID era | `Fig11-15`, `08_depth_report.md` |

### Phase 8: Innovation
| # | Script | Purpose | Output |
|---|--------|---------|--------|
| 11 | `11_innovation.py` | Counterfactual treatment recommender, AJCC staging migration, SES interaction, CONSORT diagram | `Fig16-19` (partial) |

### Phase 9: HCC vs ICC
| # | Script | Purpose | Output |
|---|--------|---------|--------|
| 12 | `12_hcc_vs_icc.py` | Full stratified analysis: baseline, KM, Cox, age×surgery type, formal interaction tests, ML by cancer type | `Fig20_HCCvsICC`, `10_hcc_vs_icc_report.md` |

### Phase 10: JAMA Surgery Level
| # | Script | Purpose | Output |
|---|--------|---------|--------|
| 13 | `13_jama_level.py` | E-value, instrumental variable, RCS splines, frailty surrogate index, geographic variation, Fine-Gray competing risk, tipping point, leave-one-out, model specification robustness | `Fig21-23`, `12_jama_level_report.md` |

### Phase 11: Model Fixes & Validation
| # | Script | Purpose | Output |
|---|--------|---------|--------|
| 14 | `14_deepsurv_fix.py` | Diagnose and fix DeepSurv loss (sort-by-risk → sort-by-time), compare architectures | `Fig24_DeepSurv_Diagnostics` |
| 15 | `15_calibration.py` | Brier scores, bootstrap calibration, risk stratification by model | `Fig25_Calibration` |
| 16 | `16_ensemble_with_deepsurv.py` | 4-model ensemble (Cox+RSF+XGB+DeepSurv) with z-score standardization, 5-fold CV | Final ensemble C-index |

---

## Key Model Results

| Model | C-index |
|-------|---------|
| **4-Model Ensemble** | **0.756 ± 0.003** |
| DeepSurv | 0.751 |
| XGBoost | 0.746 |
| Cox PH | 0.739 |
| RSF | 0.736 |

| Brier Score | Cox PH | RSF |
|-------------|--------|-----|
| 12 months | 0.165 | 0.157 |
| 36 months | 0.100 | 0.092 |
| 60 months | 0.066 | 0.059 |

---

## External Validation

| Cohort | N | C-index | 95% CI |
|--------|---|---------|--------|
| SEER Internal (5-fold CV) | 76,110 | 0.756 | — |
| TCGA-LIHC | 269 | 0.595 | 0.536–0.654 |
| ICGC-LIRI-JP | 260 | 0.547 | 0.486–0.608 |

---

## Key Findings

1. Surgery-vs-none: PSM-adjusted HR 0.35, E-value 4.3
2. Within surgical subgroup: segmental (HR 0.23) and larger (HR 0.24) — similar adjusted HRs
3. Surgery × cancer type interaction: p=0.206 (not significant)
4. Pattern holds across age ≥70, HCC, and ICC
5. ICGC validation did not succeed (CI crosses 0.5)

---

## Review History

| Stage | Reviewers | Outcome |
|-------|-----------|---------|
| Self-Review | 1 | 6 Major → all fixed |
| Panel Review | 3 (Epi + Clinical + ML) | 3 Critical + 5 Major → all fixed |
| JAMA Surgery Sim | Editor + 2 peers | Desk Reject → revised |
| Final Gates | 4 deterministic | ALL PASS |

---

## Figure Index (29 figures)

| Fig | Title |
|-----|-------|
| Fig1_KM | 6-panel Kaplan-Meier |
| Fig2_Forest | Multivariate Cox forest plot |
| Fig5_ExternalValidation | SEER→TCGA+ICGC validation |
| Fig6_CompositeAnalysis | CV + PSM + HBI + Nomogram |
| Fig7_LASSO | LASSO variable selection |
| Fig8_DCA | Decision curve analysis |
| Fig10_Nomogram | Clinical nomogram |
| Fig11_Subgroup | Subgroup interaction forest |
| Fig12_Landmark | Landmark conditional survival |
| Fig13_AFP | AFP-stratified KM |
| Fig14_Race | Race-stratified KM |
| Fig15_COVID | COVID era trends |
| Fig20_HCCvsICC | HCC vs ICC surgical comparison |
| Fig21_RCS_AgeSpline | RCS age-surgery benefit |
| Fig22_Frailty | Frailty Surrogate Index |
| Fig23_Geography | Geographic practice variation |
| Fig25_Calibration | Model calibration |
| FigS1–S4 | SHAP, Competing Risk, Age-Surgery, Temporal Trends |
