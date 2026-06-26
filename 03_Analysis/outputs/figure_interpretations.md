# Figure Interpretations — Hepatobiliary Cancer Analysis

---

## MAIN FIGURES

### Fig1_KM — Kaplan-Meier Survival Curves (6 panels)

**A. OS: Surgery vs Non-surgery** — Any surgery dramatically improves OS over no surgery (median 28-42m vs 4m). Log-rank p<0.0001.

**B. CSS: Surgery vs Non-surgery** — Similar dramatic separation for cancer-specific survival (median 14m non-surgery vs 74m surgery). Confirms that surgery primarily prevents cancer death, not just all-cause mortality.

**C. OS: By Surgery Type** — Transplant best (42m), Segmental second (30m), Local Destruction (28m), Larger Resection (24m), None (4m). **Key clinical message**: Segmental resection outperforms larger resection in median OS.

**D. OS: By AJCC Stage** — Expected stage gradient: Stage I best survival, Stage IV worst. All curves well-separated with no crossing.

**E. OS: Surgery Patients by Age** — Even among surgical patients, older age predicts worse survival. Age 80+ still derives substantial benefit (median 15m vs 2m non-surgery).

**F. CSS: HCC vs ICC vs Other** — HCC has best survival, ICC intermediate, Other histologies worst. ICC 1-year CSS ~33% vs HCC ~46%.

**Key takeaway**: Surgery benefits all subgroups; segmental resection median OS (30m) exceeds larger resection (24m).

---

### Fig2_Forest — Multivariate CSS Cox Forest Plot

All 15 covariates plotted with HR and 95% CI on log scale.

**Strongest protective factors**: Transplant (HR ~0.15), Segmental (HR ~0.23), Larger (HR ~0.24), Radiation (HR ~0.52), Chemotherapy (HR ~0.57).

**Strongest risk factors**: Stage IV (HR ~1.74), Stage III (HR ~1.67), ICC (HR ~1.29).

**Key finding**: Segmental and Larger resection have nearly identical HRs (0.23 vs 0.24) — width and position of CIs overlap almost completely.

---

### Fig5_ExternalValidation — SEER→TCGA+ICGC

**Panel A**: C-index bar chart comparing SEER internal vs TCGA vs ICGC for Cox/RSF/XGBoost. Internal C ~0.69-0.75, TCGA ~0.57-0.60, ICGC ~0.52-0.55.

**Panel B**: TCGA risk stratification — SEER-trained RSF successfully separates high/low risk groups in TCGA (p<0.05), despite lower absolute C-index.

**Panel C**: Survival curves for SEER surgery vs TCGA vs ICGC — TCGA and ICGC patients have much better baseline survival (all surgical, younger, earlier stage).

**Panel D**: Feature variance comparison — surgery_any has near-zero variance in external cohorts (the "Surgery Feature Paradox") explaining most of the C-index drop.

**Key finding**: Model generalizes to TCGA (US-based) but fails in ICGC (Japanese). Gap explained by predictor range restriction + population differences.

---

### Fig6_CompositeAnalysis — Multi-Panel Summary

**A**: 5-fold CV model comparison — Ensemble (0.756) > XGBoost (0.746) > DeepSurv (0.751) > Cox (0.739) > RSF (0.736) — all within narrow range, demonstrating stable model performance.

**B**: PSM Love Plot — SMD before matching (red dots) widely scattered; after matching (green squares) all SMD<0.1, confirming excellent covariate balance across 6,434 matched pairs.

**C**: HBI 3-tier risk stratification — Low-risk (median OS 22m), Intermediate (5m), High-risk (2m). Curves well-separated with no crossing.

**D**: Age-dependent surgery benefit — Segmental (green solid) and Larger (red dashed) HRs track closely across all ages 65-90. Both <1.0 at all ages. At age 70+, segmental curve lies on or below larger resection curve.

**Key finding**: The HBI score provides 11-fold survival difference between strata; age-surgery benefit curves confirm segmental ≈ larger at all ages.

---

### Fig7_LASSO — LASSO Variable Selection

**Left panel**: LASSO coefficient path — as penalty (alpha) increases, coefficients shrink to zero. Vertical red line marks optimal alpha via 5-fold CV.

**Right panel**: Top 10 LASSO coefficients — surgery_any dominates (largest absolute coefficient), followed by stage variables, chemotherapy, and radiation.

**Key finding**: LASSO confirms surgery as the single strongest predictor, consistent with Cox and SHAP results. Non-zero coefficients match domain knowledge.

---

### Fig8_DCA — Decision Curve Analysis

Net benefit (y-axis) vs threshold probability (x-axis, 0.01-0.50).

**RSF Model (blue)**: Positive net benefit across threshold range 0.01-0.35.

**Treat All (orange dashed)**: Highest net benefit at very low thresholds, drops steeply.

**Treat None (gray dotted)**: Reference line at zero.

**Key finding**: The RSF model provides net benefit over "treat all" or "treat none" strategies for threshold probabilities between 0.05-0.35, supporting clinical utility for risk-based treatment decisions.

---

### Fig10_Nomogram — Clinical Nomogram

Horizontal bar chart showing point contributions of 9 clinical variables to the survival prediction model.

**Largest contributors**: Surgery receipt (−10 pts, protective), Chemotherapy (−5 pts, protective), Stage IV (+4 pts, risk), Stage III (+3 pts, risk).

**Key finding**: Simple 9-variable scoring system can estimate individualized survival probability at bedside.

---

### Fig11_Subgroup — Subgroup Interaction Forest Plot

14 subgroup estimates of surgery CSS HR with 95% CIs.

All point estimates < 1.0 (surgery protective in every subgroup). Tightest CIs in HCC, Stage I-II. Wider CIs in ICC, Stage IV (smaller sample sizes). No subgroup shows qualitative interaction (all favor surgery).

**Key finding**: Surgery benefit is universal — no subgroup where surgery is harmful. Age 80+ still shows HR ~0.30.

---

### Fig12_Landmark — Landmark Conditional Survival

Three panels showing survival from landmark time points (0, 12, 24 months), conditioned on having survived to that point.

At each landmark, surgery patients (blue) maintain higher survival than non-surgery (orange). The gap narrows but persists: at 24-month landmark, surgical patients still have ~2× survival probability.

**Key finding**: Surgery benefit is durable — even after surviving 2 years, surgical patients continue to have better prognosis.

---

### Fig13_AFP — AFP- Stratified Survival

Two panels: AFP Negative (left) vs Positive (right).

Surgery benefit present in both groups. AFP+ patients have worse baseline survival (non-surgery median ~3m) but still derive proportional benefit from surgery.

**Key finding**: AFP positivity is prognostic (worse survival) but not predictive (surgery benefit similar magnitude in both groups).

---

### Fig14_Race — Race-Stratified Survival

Four panels: NHW, NHB, NHAPI, Hispanic.

Surgery benefit observed across all racial groups. NHAPI shows best survival (both surgical and non-surgical) — consistent with known "Asian paradox" in HCC. NHB shows wider CIs due to smaller sample.

**Key finding**: Surgery benefits all racial groups. NHAPI survival advantage persists after adjustment, likely reflecting different HCC etiology (HBV-dominant vs HCV/alcohol).

---

### Fig15_COVID — COVID Era Impact

Three panels: Annual case counts, Surgery rate, Median OS.

Case counts dipped slightly in 2020, recovered by 2021-2022. Surgery rate dropped from ~24% (2019) to ~21% (2020). Median OS dropped from 12m (2019) to 6m (2020), with partial recovery to 10m (2021). Red shading = COVID period.

**Key finding**: COVID-19 pandemic temporarily reduced surgery rates and OS, with incomplete recovery by 2022.

---

## SUPPLEMENTARY FIGURES

### FigS1_SHAP — SHAP Feature Importance

Beeswarm plot showing SHAP values for top 15 features from XGBoost model.

surgery_any dominates (largest spread of SHAP values), followed by age_c, stage_4, stage_3, chemotherapy. Red points = high feature value, blue = low. For surgery_any, red points cluster on left (protective) — consistent with clinical expectation.

**Key finding**: SHAP confirms surgery as the single most influential predictor, with direction and magnitude consistent with Cox and LASSO results.

---

### FigS2_CompetingRisk — Competing Risk CIF

Three panels: Cumulative Incidence Functions by (A) surgery status, (B) age band, (C) AJCC stage.

Solid lines = cancer-specific death CIF; dashed = other-cause death CIF. Cancer death dominates (~82% of all deaths). Other-cause death CIF higher in older age bands (expected).

**Key finding**: Cancer-specific death is the dominant competing event; other-cause death contributes modestly even in the oldest age group.

---

### FigS3_AgeSurgeryBenefit — Age-Dependent Surgery Benefit Spline

Surgery CSS HR plotted against continuous age (65-90). Blue dots = unadjusted, Red line = adjusted spline.

HR is lowest (~0.22) at age ~72, gradually rises to ~0.35 at age 90. Benefit inflection accelerates after age 86. All HRs remain <0.40.

**Key finding**: Surgery benefit attenuates with age but remains clinically significant at all ages (HR <0.40 at age 90). Inflection at age 86 where attenuation accelerates.

---

### FigS4_TemporalTrends — Temporal Trends (2004-2022)

Dual Y-axis: Surgery rate (blue, left axis) and Median OS (red, right axis).

Surgery rate: ~14% (2004) → ~24% (2019) → ~21% (2020-2022). Median OS: 3m (2004) → 12m (2019) → 6-10m (2020-2022). COVID dip visible in both series.

**Key finding**: Steady improvement in both surgery utilization and median OS over 18 years, disrupted by COVID-19 pandemic with incomplete recovery.

---

### Fig20_HCCvsICC — HCC vs ICC Surgical Strategy Comparison

Four panels directly comparing HCC and ICC.

**A-B**: Age-dependent surgery HR for Segmental vs Larger in HCC and ICC separately — both show segmental HR ≤ larger HR at ages ≥70.

**C**: Surgery type distribution — HCC has higher transplant rate, ICC has higher local destruction rate.

**D-E**: Stage-specific surgery HR — transplant is dramatically better in HCC than ICC (HR 0.12 vs 0.24). Segmental and Larger are similar in both.

**Key finding**: The segmental≈larger pattern observed in HCC extends to ICC. Transplant benefit is cancer-type dependent (stronger in HCC).

---

### Fig21_RCS_AgeSpline — Restricted Cubic Spline Age Effect

Smooth spline curve showing non-linear relationship between age and surgery HR. Peak benefit at age ~72. Plateau from 65-82. Gradual attenuation from 82-90.

**Key finding**: Surgery benefit is not linear with age — it has an optimal window (65-82) and accelerates in decline after 86.

---

### Fig22_Frailty — Frailty Surrogate Index

Three panels: KM curves by FSI group (Fit / Pre-frail / Frail). Within each FSI group, surgery vs non-surgery comparison.

Even in the Frail group (age ≥80 + stage ≥III), surgery provides absolute 13-month median OS benefit (15m vs 2m). Surgery rate decreases with frailty (31% → 15% → 10%).

**Key finding**: Frailty should not automatically disqualify patients from surgery — even frail patients derive substantial survival benefit.

---

### Fig23_Geography — Geographic Practice Variation

Bar chart comparing surgery rate, chemo rate, and median OS across 4 SES quartile regions.

Higher-income regions have modestly higher surgery rates (24% vs 20%) and median OS (8m vs 6m). Segmental/Larger ratio is consistently 3-4× across all regions.

**Key finding**: Parenchymal-sparing practice pattern is consistent nationally; SES gradient exists in access to surgery but is modest in magnitude.

---

### Fig25_Calibration — Model Calibration & Validation

Six panels:

**A-C**: Risk stratification by Cox (A), RSF (B), XGBoost (C). All three models cleanly separate Low/Medium/High risk groups with well-ordered survival curves.

**D**: Calibration plot — observed vs predicted event probability at 36 months. Points cluster near diagonal, confirming good calibration.

**E**: C-index comparison — XGBoost (0.746) > RSF (0.736) > Cox (0.739). All within narrow range.

**F**: Brier score over time — RSF consistently lower (better) than Cox PH at all time points.

**Key finding**: Models are well-calibrated (not just good discrimination). RSF has best calibration among tested models.
