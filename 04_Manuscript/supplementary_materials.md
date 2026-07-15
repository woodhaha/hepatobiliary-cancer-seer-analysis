# Supplementary Materials

## eFigure 1: CONSORT Flow Diagram

```
SEER Research Plus Data (Nov 2024 submission)
  N = 171,286 hepatobiliary cancer records
         │
         ├── Exclude: Age < 18 (n = 1,399)
         │
         ▼
  N = 169,887 adult patients
         │
         ├── Exclude: Non-hepatobiliary site (n = 0)
         │   └── Liver (C22.0) + Intrahepatic Bile Duct (C22.1)
         │
         ▼
  N = 169,887
         │
         ├── Exclude: Missing survival data (n = 5,504)
         │
         ▼
  N = 164,383 with complete survival
         │
         ├── Study Cohort: Age ≥ 65 years
         │   N = 76,110
         │      │
         │      ├── Non-surgery: n = 59,821 (78.6%)
         │      │   └── Included in: Baseline, Cox, PSM
         │      │
         │      └── Surgery: n = 16,289 (21.4%)
         │           ├── Local Destruction: n = 10,771 (14.2%)
         │           ├── Segmental Resection: n = 2,486 (3.3%)
         │           ├── Larger Resection: n = 699 (0.9%)
         │           ├── Transplant: n = 1,778 (2.3%)
         │           └── Other/NOS: n = 555 (0.7%)
         │
         ├── PSM Analysis
         │   └── 6,434 matched pairs (SMD 0.190 → 0.030)
         │
         ├── ML Training/Testing
         │   ├── Train (2004-2017): n = 53,277
         │   └── Test (2018-2022): n = 22,833
         │
         └── External Validation
             ├── TCGA-LIHC: n = 269 (all surgical HCC, US)
             └── ICGC-LIRI-JP: n = 260 (all surgical HCC, Japan)
```

---

## eFigure 2: PSM Love Plot

**Before Matching**: Mean SMD = 0.190
- Stage IV: SMD = 0.72
- Chemotherapy: SMD = 0.36
- Age: SMD = 0.33
- Stage I: SMD = 0.31
- Radiation: SMD = 0.22

**After Matching** (1:1 nearest-neighbor, caliper=0.05): Mean SMD = 0.030
- All covariates SMD < 0.10 (adequately balanced)
- Largest residual: Income SMD = 0.06

Full love plot: `03_Analysis/figures/Fig6_CompositeAnalysis.png` (Panel B)

---

## eFigure 3: 5-Fold CV Model Comparison

| Model | Fold 1 | Fold 2 | Fold 3 | Fold 4 | Fold 5 | Mean ± SD |
|-------|--------|--------|--------|--------|--------|-----------|
| Cox PH | 0.740 | 0.738 | 0.739 | 0.740 | 0.738 | 0.739 ± 0.001 |
| RSF | 0.735 | 0.737 | 0.738 | 0.734 | 0.736 | 0.736 ± 0.002 |
| XGBoost | 0.745 | 0.748 | 0.746 | 0.744 | 0.747 | 0.746 ± 0.002 |
| **Ensemble** | 0.756 | 0.754 | 0.757 | 0.755 | 0.758 | 0.756 ± 0.003 |

Full figure: `03_Analysis/figures/Fig6_CompositeAnalysis.png` (Panel A)

---

## eFigure 4: SHAP Feature Importance

Top 10 features by SHAP value (XGBoost model):

| Rank | Feature | Mean SHAP |
|------|---------|-----------|
| 1 | surgery_any | 0.128 |
| 2 | stage_4 | 0.065 |
| 3 | stage_3 | 0.040 |
| 4 | age_c | 0.035 |
| 5 | chemotherapy | 0.031 |
| 6 | is_icc | 0.025 |
| 7 | stage_2 | 0.018 |
| 8 | radiation | 0.015 |
| 9 | cirrhosis | 0.012 |
| 10 | grade_poor | 0.009 |

Full figure: `03_Analysis/figures/FigS1_SHAP.png`

---

## eFigure 5: External Validation (TCGA-LIHC + ICGC-LIRI-JP)

| Cohort | N | Surgery% | Age | Dead% | Med OS | Cox C | RSF C | XGB C |
|--------|---|----------|-----|-------|--------|-------|-------|-------|
| SEER Train | 53,277 | 21% | 74 | 83% | 7m | 0.739 | 0.736 | 0.746 |
| SEER Test | 22,833 | 22% | 74 | 83% | 8m | — | — | — |
| TCGA-LIHC | 269 | 100% | 59 | 7% | 22m | 0.595 | 0.567 | 0.592 |
| ICGC-LIRI-JP | 260 | 100% | 67 | 32% | 21m | 0.522 | 0.551 | 0.547 |

**Predictor Range Restriction**: In all-surgical external cohorts, `surgery_any` variance = 0 (vs SEER σ²=0.17). This eliminates the model's strongest predictor, explaining the external C-index gap.

Full figure: `03_Analysis/figures/Fig5_ExternalValidation.png`

---

## eFigure 6: Temporal Trends (2004-2022)

| Period | Surgical Rate | Transplant Rate | Median OS | Chemotherapy |
|--------|--------------|-----------------|-----------|-------------|
| 2004-2009 | 21% | 2.3% | 5m | 29% |
| 2010-2015 | 22% | 2.4% | 8m | 31% |
| 2016-2019 | 23% | 2.4% | 11m | 32% |
| 2020-2022 (COVID) | 21% | 2.4% | 6m | 27% |

Full figure: `03_Analysis/figures/FigS4_TemporalTrends.png` + `Fig15_COVID.png`

---

## eFigure 7: Frailty Surrogate Index

| FSI Group | N | Age | Surgery% | Med OS (Surg) | Med OS (Non) | HR |
|-----------|----|-----|----------|---------------|--------------|-----|
| Fit | 34,248 | 71.2 | 31% | 33m | 8m | 0.27 |
| Pre-frail | 30,324 | 76.3 | 15% | 22m | 4m | 0.32 |
| Frail | 11,538 | 82.1 | 10% | 15m | 2m | 0.33 |

FSI components: age ≥80 (2pts), age ≥75 (1pt), stage ≥III (2pts), poor grade (1pt).

Full figure: `03_Analysis/figures/Fig22_Frailty.png`

---

## eFigure 8: Restricted Cubic Spline — Age-Surgery Benefit

Surgery benefit follows a U-shaped relationship with age:
- Peak benefit at age ~72 (HR ≈ 0.24)
- Stable plateau from ages 65-82
- Inflection point at age 86 — beyond this, attenuation accelerates
- At age 90, HR ≈ 0.40 (benefit reduced but still significant)

Full figure: `03_Analysis/figures/Fig21_RCS_AgeSpline.png` + `FigS3_AgeSurgeryBenefit.png`

---

## eFigure 9: Model Calibration — Brier Scores and Calibration Plots

Model calibration was assessed via time-dependent Brier scores at 12, 36, and 60 months across all four models:

| Model | 12-month | 36-month | 60-month |
|-------|----------|----------|----------|
| RSF | 0.157 | 0.092 | 0.059 |
| Cox PH | 0.165 | 0.100 | 0.066 |
| XGBoost | 0.163 | 0.097 | 0.063 |
| DeepSurv | 0.170 | 0.104 | 0.069 |

Calibration plots (eFigure 9 Panel D) confirmed good agreement between predicted and observed survival at 36 months for all models, with calibration slopes ranging from 0.91 (DeepSurv) to 0.97 (RSF).

Full figure: `03_Analysis/figures/Fig25_Calibration.png`

---

## eTable 1: Full Baseline Characteristics

| Variable | Non-surgery (N=59,821) | Surgery (N=16,289) | Total (N=76,110) | P-value |
|----------|----------------------|--------------------|-------------------|---------|
| **Age (mean±SD)** | 74.6 ± 7.0 | 72.4 ± 5.8 | 74.2 ± 6.8 | <0.001 |
| Age 65-69 | 20,747 (34.7%) | 7,625 (46.8%) | 28,372 (37.3%) | |
| Age 70-74 | 14,310 (23.9%) | 4,066 (25.0%) | 18,376 (24.1%) | |
| Age 75-79 | 11,449 (19.1%) | 2,782 (17.1%) | 14,231 (18.7%) | |
| Age 80+ | 13,315 (22.3%) | 1,816 (11.1%) | 15,131 (19.9%) | |
| **Male** | 39,012 (65.2%) | 10,721 (65.8%) | 49,733 (65.4%) | 0.155 |
| **Married** | 31,475 (52.6%) | 9,932 (61.0%) | 41,407 (54.4%) | <0.001 |
| **Race** | | | | <0.001 |
| NHW | 31,883 (53.3%) | 7,983 (49.0%) | 39,866 (52.4%) | |
| NHB | 5,754 (9.6%) | 1,028 (6.3%) | 6,782 (8.9%) | |
| NHAPI | 8,568 (14.3%) | 3,916 (24.0%) | 12,484 (16.4%) | |
| Hispanic | 11,241 (18.8%) | 2,489 (15.3%) | 13,730 (18.0%) | |
| **AJCC Stage** | | | | <0.001 |
| Stage I | 25,002 (41.8%) | 9,852 (60.5%) | 34,854 (45.8%) | |
| Stage II | 9,625 (16.1%) | 3,191 (19.6%) | 12,816 (16.8%) | |
| Stage III | 9,644 (16.1%) | 1,674 (10.3%) | 11,318 (14.9%) | |
| Stage IV | 15,550 (26.0%) | 1,572 (9.7%) | 17,122 (22.5%) | |
| **Poor Grade** | 19,558 (32.7%) | 4,464 (27.4%) | 24,022 (31.6%) | <0.001 |
| **ICC Histology** | 9,553 (16.0%) | 2,196 (13.5%) | 11,749 (15.4%) | <0.001 |
| **Chemotherapy** | 19,912 (33.3%) | 4,381 (26.9%) | 24,293 (31.9%) | <0.001 |
| **Radiation** | 8,616 (14.4%) | 1,076 (6.6%) | 9,692 (12.7%) | <0.001 |
| **Cirrhosis** | 5,035 (8.4%) | 2,095 (12.9%) | 7,130 (9.4%) | <0.001 |
| **Tumor Size (mm)** | 138 ± 261 | 155 ± 277 | 143 ± 266 | <0.001 |
| **AFP (when available)** | | | | |
| Positive | 17,523 (41.4%) | 4,481 (33.6%) | 22,004 (28.9%)* | |
| Unknown | 17,467 (29.2%) | 2,950 (18.1%) | 20,417 (26.8%)* | |
| **Income ($10K)** | 8.1 ± 2.5 | 8.5 ± 2.4 | 8.2 ± 2.5 | <0.001 |

*Percent of total cohort (including 'Unknown' category)

---

## eTable 2: Complete Multivariate Cox Regression (CSS)

| Variable | HR | 95% CI (lower) | 95% CI (upper) | P-value |
|----------|-----|----------------|----------------|---------|
| age_c | 1.02 | 1.01 | 1.02 | <0.001 |
| male | 1.05 | 1.03 | 1.07 | <0.001 |
| married | 0.97 | 0.95 | 0.99 | <0.001 |
| race_nhb | 0.91 | 0.88 | 0.94 | <0.001 |
| race_nhapi | 0.79 | 0.77 | 0.81 | <0.001 |
| race_hispanic | 0.92 | 0.90 | 0.94 | <0.001 |
| stage_2 | 1.22 | 1.19 | 1.25 | <0.001 |
| stage_3 | 1.67 | 1.62 | 1.71 | <0.001 |
| stage_4 | 1.74 | 1.70 | 1.78 | <0.001 |
| grade_poor | 1.04 | 1.02 | 1.06 | <0.001 |
| is_icc | 1.29 | 1.26 | 1.32 | <0.001 |
| surg_local_destruction | 0.26 | 0.25 | 0.27 | <0.001 |
| surg_segmental_resection | 0.23 | 0.22 | 0.25 | <0.001 |
| surg_larger_resection | 0.24 | 0.22 | 0.27 | <0.001 |
| surg_transplant | 0.15 | 0.14 | 0.17 | <0.001 |
| chemotherapy | 0.57 | 0.56 | 0.59 | <0.001 |
| radiation | 0.52 | 0.50 | 0.53 | <0.001 |
| cirrhosis | 0.82 | 0.80 | 0.85 | <0.001 |

C-index = 0.739 | Log-likelihood ratio test: p < 0.0001

---

## eTable 3: Stratified E-values by Age Band

| Age Band | Surgery HR | E-value (point) | E-value (CI-bound) | Robustness |
|----------|-----------|-----------------|-------------------|------------|
| 65-69 | 0.30 | 6.2 | 6.4 | Very Robust |
| 70-74 | 0.30 | 6.2 | 6.4 | Very Robust |
| 75-79 | 0.32 | 5.8 | 6.0 | Very Robust |
| 80+ | 0.31 | 6.0 | 6.2 | Very Robust |
| HCC | 0.25 | 7.5 | 7.7 | Very Robust |
| ICC | 0.28 | 6.6 | 6.9 | Very Robust |
| Stage I | 0.25 | 7.5 | 7.7 | Very Robust |
| Stage IV | 0.32 | 5.8 | 6.0 | Very Robust |
| PSM-adjusted | 0.35 | 5.2 | 5.5 | Very Robust |

All E-values > 5, indicating extreme robustness to unmeasured confounding.

---

## eTable 4: Model Specification Robustness

| Model | Features | Surgery HR | C-index |
|-------|----------|-----------|---------|
| Unadjusted | surgery_any | 0.27 | 0.616 |
| + Demographics | + age, sex, marital | 0.29 | 0.650 |
| + Stage | + AJCC stage | 0.31 | 0.683 |
| + Tumor | + ICC, grade | 0.31 | 0.684 |
| + Treatment | + chemo, radiation | 0.25 | 0.736 |
| Full (all 15) | + race, cirrhosis, income | 0.25 | 0.738 |

---

## eTable 5: Leave-One-Out Sensitivity Analysis

| Excluded Subgroup | Remaining N | Surgery HR | Δ from base |
|-------------------|-------------|-----------|-------------|
| None (full cohort) | 76,110 | 0.25 | — |
| Stage IV patients | 58,988 | 0.29 | +0.04 |
| Age 80+ patients | 60,979 | 0.25 | 0.00 |
| ICC patients | 64,361 | 0.25 | 0.00 |
| Age 80+ | 58,582 | 0.30 | +0.05 |
| Low income (<median) | 42,299 | 0.30 | +0.05 |
| Early era (2004-2010) | 51,410 | 0.27 | +0.02 |

---

## eTable 6: Geographic Practice Variation

| SES Quartile | N | Surgery% | Transplant% | Chemo% | Median OS | Segmental/Larger Ratio |
|-------------|----|----------|------------|--------|-----------|----------------------|
| Q1 (Lowest) | 23,580 | 20% | 2.4% | 31.8% | 6m | 3.4 |
| Q2 | 17,008 | 20% | 2.2% | 32.0% | 7m | 3.6 |
| Q3 | 16,949 | 22% | 2.2% | 30.5% | 7m | 3.2 |
| Q4 (Highest) | 18,573 | 24% | 2.5% | 33.2% | 8m | 4.1 |

---

---

## eTable 7: Model Hyperparameters

| Parameter | Cox PH | RSF | XGBoost | DeepSurv |
|-----------|--------|-----|---------|----------|
| **Implementation** | sksurv 0.22 | sksurv 0.22 | xgboost 2.1 | PyTorch 2.5 |
| **Architecture** | Linear | 100 trees | Gradient boosting | 17→64→32→16→1 |
| **Regularization** | L2 (alpha=0.01) | — | λ=1.0 | weight_decay=1e-5 |
| **Min samples split** | — | 50 | — | — |
| **Min samples leaf** | — | 20 | — | — |
| **Max depth** | — | None | 4 | — |
| **Max features** | — | sqrt | — | — |
| **Learning rate** | — | — | 0.05 | 1e-3 |
| **Subsample** | — | — | 0.8 | — |
| **Colsample bytree** | — | — | 0.8 | — |
| **Min child weight** | — | — | 5 | — |
| **Boosting rounds** | — | — | 150 | — |
| **Batch size** | — | — | — | 256 |
| **Epochs (max)** | — | — | — | 200 |
| **Early stopping** | — | — | — | patience=15 |
| **Optimizer** | — | — | — | Adam |
| **Dropout** | — | — | — | 0.3 |
| **BatchNorm** | — | — | — | Yes |
| **Loss function** | Cox partial likelihood | Log-rank splitting | Cox partial likelihood (Efron) | Cox partial likelihood |
| **Random seed** | — | 42 | 42 | — |

---

## eTable 8: Schoenfeld Residual Tests

| Variable | Test Statistic | P-value | PH Violation |
|----------|---------------|---------|:---:|
| age_c | 89.2 | <0.001 | ⚠ |
| male | 0.02 | 0.881 | ✓ |
| married | 6.1 | 0.014 | ⚠ |
| stage_2 | 28.4 | <0.001 | ⚠ |
| stage_3 | 45.1 | <0.001 | ⚠ |
| stage_4 | 15.3 | <0.001 | ⚠ |
| grade_poor | 2.7 | 0.099 | ✓ |
| is_icc | 9.9 | 0.002 | ⚠ |
| chemotherapy | 112.7 | <0.001 | ⚠ |
| radiation | 67.3 | <0.001 | ⚠ |
| cirrhosis | 1.5 | 0.216 | ✓ |
| surgery_any | 58.6 | <0.001 | ⚠ |

*Note: With N=76,110, formal PH tests detect trivially small violations due to high statistical power. Clinically, surgery exhibits time-varying effects (early perioperative risk followed by long-term survival benefit). Landmark analyses at 12 and 24 months corroborate the time-averaged HRs reported in the primary analysis.*

---

## Supplementary References

S1. TCGA Research Network. Comprehensive and integrative genomic characterization of hepatocellular carcinoma. *Cell*. 2017;169(7):1327-1341.
S2. Totoki Y, et al. Trans-ancestry mutational landscape of hepatocellular carcinoma genomes. *Nat Genet*. 2014;46:1267-1273.
S3. VanderWeele TJ, Ding P. Sensitivity analysis in observational research: introducing the E-value. *Ann Intern Med*. 2017;167(4):268-274.
S4. Baiocchi M, Cheng J, Small DS. Instrumental variable methods for causal inference. *Stat Med*. 2014;33(13):2297-2340.

---

## eFigure 10: Composite Analysis — PSM + HBI + Landmark + Age-Surgery

**Panel A**: 5-Fold Cross-Validated Model Comparison (C-index by fold for Cox, RSF, XGBoost, DeepSurv, Ensemble)
**Panel B**: PSM Love Plot (SMD before/after matching, 0.190→0.030)
**Panel C**: HBI Risk Score Stratification (Low/Intermediate/High)
**Panel D**: Age-Surgery Benefit by FSI Group

Full figure: `03_Analysis/figures/Fig6_CompositeAnalysis.png`

---

## eFigure 11: Clinical Nomogram for Individualized Survival Prediction

A 7-variable nomogram predicting 36-month CSS probability based on: surgery type, AJCC stage, age, ICC histology, chemotherapy, radiation, and cirrhosis.

- C-index: 0.739 (Cox PH)
- Calibration-in-the-large: good agreement at 36 months

Full figure: `03_Analysis/figures/Fig10_Nomogram.png`

---

## eTable 9: Clinical Decision Matrix

| Age | Cancer Type | Stage | FSI Group | Recommended Surgery | Rationale |
|-----|-------------|-------|-----------|-------------------|-----------|
| 65-69 | HCC | I-II | Fit | Segmental or Larger | Similar outcomes; larger if margin concern |
| 65-69 | HCC | III-IV | Fit | Segmental | Palliative intent; minimize morbidity |
| 65-69 | ICC | I-II | Fit | Segmental | Segmental HR 0.22 vs Larger HR 0.26 |
| 70-74 | HCC/ICC | Any | Fit/Pre-frail | Segmental | Equivalent survival; lower morbidity |
| 75-79 | HCC/ICC | Any | Fit/Pre-frail | Segmental | Segmental superior in all age bands ≥75 |
| 80+ | HCC/ICC | Any | Pre-frail/Frail | Segmental | Benefit persists; 13-month OS gain |
| Any | Any | IV | Any | Non-surgical | Systemic therapy first; palliation |
| Any | Any | I-II | Frail | Segmental/Local | Frail patients gain 13 months OS |
