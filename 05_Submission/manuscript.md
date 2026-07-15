# Surgical Resection in Elderly Hepatobiliary Cancer: A Population-Based Cohort Study Comparing Parenchymal-Sparing and Major Hepatectomy

**Qigang Xu**$^{2,\dagger}$ · **Yongyu Bai**$^{1,\dagger}$ · **Yiqi Cai**$^{1}$ · **Zhuha Zhou**$^{1,*}$

$^{1}$ Department of Gastroenterology Surgery, The First Affiliated Hospital of Wenzhou Medical University, Zhejiang, China
$^{2}$ Department of Hepatobiliary and Pancreatic Surgery, The First Affiliated Hospital of Wenzhou Medical University, Zhejiang, China
$\dagger$ Equal contribution · $^{*}$ Corresponding author: zhouzhuha@wmu.edu.cn

**Running head**: Segmental vs Major Hepatectomy in Elderly Hepatobiliary Cancer

**Key words**: hepatocellular carcinoma, intrahepatic cholangiocarcinoma, elderly, SEER, segmental resection, E-value, external validation

**Word count**: Abstract 249 | Text ~2,600 | Figures 6 | Tables 3 | References 25

---

## ABSTRACT

**Importance**: The survival benefit of surgical resection in elderly hepatobiliary cancer is well-established, but whether segmental/wedge resection offers comparable outcomes to major hepatectomy across both HCC and intrahepatic cholangiocarcinoma (ICC) is unknown.

**Objective**: To quantify the survival benefit of surgical resection and compare outcomes by extent of resection.

**Design, Setting, and Participants**: Population-based cohort study using SEER 18-registry data (2004–2022), with external validation on TCGA-LIHC and ICGC-LIRI-JP. Participants were 76,110 patients aged ≥65 with histologically confirmed hepatobiliary malignancy.

**Exposures**: Five surgical categories: none, local destruction, segmental/wedge resection, larger resection (lobectomy/extended lobectomy), and liver transplantation.

**Main Outcomes and Measures**: Cancer-specific survival. Multivariate Cox regression with propensity score matching, E-value analysis, and a 4-model machine learning ensemble.

**Results**: Surgical resection was strongly protective versus no surgery (PSM-adjusted HR 0.35; 95% CI 0.33–0.37; E-value 5.2). Within the surgical subgroup, segmental resection (HR 0.23; 95% CI 0.22–0.25) and larger resection (HR 0.24; 95% CI 0.22–0.27) showed similar adjusted survival. However, this within-surgery comparison is observational and susceptible to confounding (E-value ~1.3); liver function, tumor anatomy, and performance status — absent from SEER — could plausibly account for the observed similarity. The surgery × cancer type interaction was not significant (p=0.206). Median OS ranged from 4 months (non-surgery) to 42 months (transplant). The 4-model ensemble achieved C-index 0.756 internally. External validation was limited: TCGA-LIHC C-index 0.595 (95% CI 0.536–0.654); ICGC-LIRI-JP C-index 0.547 (95% CI 0.486–0.608), essentially chance-level.

**Conclusions and Relevance**: Surgical resection confers a robust survival benefit in elderly hepatobiliary cancer patients. The comparison between segmental and larger resection is observational and sensitive to unmeasured confounding; these findings should be interpreted as hypothesis-generating rather than evidence of equivalence.

---

## INTRODUCTION

Hepatobiliary cancers, comprising hepatocellular carcinoma (HCC) and intrahepatic cholangiocarcinoma (ICC), represent some of the most lethal malignancies worldwide [1]. With global population aging, the proportion of elderly patients (age ≥65 years) diagnosed with these cancers continues to rise annually [2]. For elderly patients with resectable hepatobiliary cancer, multiple surgical options exist—from segmental or wedge resection to lobectomy, extended lobectomy, and liver transplantation [3]. However, controversy persists: more radical surgeries may offer superior oncologic clearance but carry higher perioperative risk, particularly in elderly patients with comorbidities [4]. Current guidelines [18,19] recommend surgical resection based on tumor stage and liver function, with BCLC staging providing the framework for treatment allocation [20], yet evidence specifically guiding the choice between segmental and major resection in elderly patients remains limited.

Zhang et al. first reported from SEER data that for 2,371 elderly patients with stage I–II HCC, segmental or wedge resection provided survival outcomes comparable to larger resection when transplantation was unavailable [5]. However, this analysis was limited to early-stage HCC (2004–2011) and did not include ICC. Xia et al. subsequently compared surgical types in 10,174 HCC patients but combined wedge and segmental resection as a single category without age-specific stratification of surgical subtypes [6]. For ICC, an XGBoost model predicted survival in 1,055 surgical patients but did not compare surgical extents [7]. No prior study has simultaneously compared segmental versus larger resection outcomes in elderly patients with both HCC and ICC using population-based data with formal causal sensitivity analyses and external validation.

To address these gaps, we conducted a comprehensive analysis with three objectives: (1) to quantify the survival benefit of surgical resection in elderly hepatobiliary cancer using robust causal sensitivity analyses; (2) to compare adjusted outcomes by surgical extent within the operated subgroup across HCC and ICC; and (3) to assess robustness through E-value analysis, propensity score matching, and multiple sensitivity analyses.

---

## METHODS

### Data Source and Study Population

The SEER Research Plus Data (November 2024 submission) was queried for patients with hepatobiliary cancer (primary site: liver C22.0 or intrahepatic bile duct C22.1) diagnosed between January 1, 2004 and December 31, 2022. SEER covers approximately 28% of the US population across 18 registries [8]. Inclusion criteria were: (1) age ≥65 years at diagnosis; (2) histologically confirmed hepatobiliary malignancy (ICD-O-3: 8170–8175, 8180 for HCC; 8160–8162 for ICC); (3) known AJCC stage and survival data; and (4) complete data on key covariates (age, sex, marital status, race, stage, grade, cancer type, surgery, chemotherapy, radiation, cirrhosis, tumor size). Of 171,286 records, 76,110 (44.4%) met inclusion criteria; exclusions included age <18 (n=1,399), missing survival data (n=5,504), age <65 (n=87,773), and incomplete covariate data (n=500). A comparison of included versus excluded patients on available demographics is provided (eTable 10). Patients with other hepatobiliary histologies (n=6,981; 9.2%) including combined HCC-ICC, hepatoblastoma, and sarcoma NOS were retained in the full cohort for primary analyses; a sensitivity analysis excluding these patients is reported (eTable 10). The study was exempt from institutional review board approval as SEER, TCGA, and ICGC data are publicly available and de-identified.

### Statistical Power

With N=76,110 and approximately 82.7% mortality, the study has >99% power at α=0.05 to detect a hazard ratio of 1.1 or smaller for the surgery-versus-none comparison. For the surgical subgroup comparison (segmental vs larger resection, n≈3,185), the study has 80% power to detect a HR of 1.25 or greater. This is adequate for the primary comparison; subgroup analyses (e.g., within individual age-cancer type strata) have proportionally lower power and are interpreted with caution.

### Variable Definitions

The primary exposure was surgical treatment categorized into five groups: (1) None (no cancer-directed surgery); (2) Local Destruction (radiofrequency ablation, percutaneous ethanol injection, cryoablation; SEER codes 10–27); (3) Segmental Resection (wedge or segmental resection; codes 30–38); (4) Larger Resection (lobectomy, extended lobectomy; codes 50–59); and (5) Liver Transplantation (codes 60–66).

Due to changes in AJCC staging across the study period (3rd, 6th, 7th editions, and EOD 2018), a unified staging algorithm was implemented that harmonized stage groups to I–IV by prioritizing the most granular staging system for each diagnosis year, consistent with prior SEER-based hepatobiliary cancer studies [5,6].

Primary outcomes were overall survival (OS) and cancer-specific survival (CSS).

### Statistical Analysis

Baseline characteristics were compared using Student's t-test (continuous) and χ² test (categorical). Kaplan-Meier curves were generated with log-rank tests. Univariate and multivariate Cox proportional hazards regression were conducted for CSS. Stratified analyses were performed by age band (65–69, 70–74, 75–79, 80+), AJCC stage (I–IV), and cancer type (HCC, ICC). Formal interaction testing was conducted between surgery and cancer type, age, and stage.

### Propensity Score Matching

1:1 propensity score matching was performed using logistic regression with nearest-neighbor matching (caliper=0.05) [9]. Covariates included age, sex, marital status, stage, grade, cancer type, chemotherapy, radiation, and tumor size. Standardized mean differences (SMD) were calculated before and after matching; SMD <0.1 was considered adequate balance [10].

### Machine Learning Models

Four survival models were trained: L2-regularized Cox PH, Random Survival Forest (100 trees, min_samples_leaf=20) [11], XGBoost (survival:cox objective, max_depth=4, learning_rate=0.05) [11], and DeepSurv (3-layer MLP: 17→64→32→16→1 with ReLU activations, batch normalization, dropout=0.3, trained with the negative log partial likelihood of the Cox model; Adam optimizer, lr=1e-3, weight_decay=1e-5, batch_size=256, 200 epochs with early stopping patience=15). Models were trained using temporal split (2004–2017 train, 2018–2022 test). A 5-fold cross-validated 4-model ensemble combined Cox, RSF, XGBoost, and DeepSurv predictions via z-score standardization and equal-weight averaging of out-of-fold risk scores. Performance was assessed using Harrell's C-index with Uno's concordance correction [12]. Feature importance was evaluated using permutation importance and SHAP values [10]. Model calibration was assessed at 12, 36, and 60 months using time-dependent Brier scores and calibration-in-the-large [13].

### External Validation

Two independent cohorts served for external validation: TCGA-LIHC (n=269, all-surgical HCC, US multi-center) [14] and ICGC-LIRI-JP (n=260, all-surgical HCC, Japanese) [15]. SEER-trained models were applied to external cohorts without retraining.

### Sensitivity Analyses

(1) E-value analysis for unmeasured confounding [16]; (2) Cause-specific hazard regression and subdistribution hazard (Fine-Gray) models for competing risks [17]; (3) Leave-one-out analysis excluding key subgroups, including exclusion of non-HCC/ICC histologies; (4) Sequential model specification testing; (5) Restricted cubic splines for non-linear age effects; (6) Geriatric Risk Index (GRI) using age, stage, and grade; (7) Landmark conditional survival analysis at 12 and 24 months; (8) Schoenfeld residual tests for proportional hazards assumption with time-varying coefficient estimates.

### Temporal, Geographic, and Clinical Tool Development

Annual trends (2004–2022) were analyzed. Geographic variation was assessed across socioeconomic quartiles. COVID-19 era (2020–2022) was compared to the pre-COVID period (2017–2019). A 7-variable Hepatobiliary Cancer Index (HBI) risk score was derived and validated, with a clinical nomogram and decision matrix constructed.

This study follows the TRIPOD (Transparent Reporting of a multivariable prediction model for Individual Prognosis Or Diagnosis) statement for prediction model development with external validation (Type 3). All statistical tests were two-sided with α=0.05. Analyses used Python 3.13 (scikit-survival v0.22, lifelines v0.28, xgboost v2.1, PyTorch 2.5) and R 4.4.

---

## RESULTS

### Study Population

**eFigure 1** summarizes cohort selection. The analytic cohort comprised 76,110 elderly patients: 57,380 HCC (75.4%), 11,749 ICC (15.4%), and 6,981 other hepatobiliary malignancies (9.2%; including combined HCC-ICC, hepatoblastoma, and sarcoma NOS). Median age was 74 years (IQR 68–80); 65.4% were male. Median follow-up was 49 months (IQR 16–96) for OS and 74 months (IQR 24–120) for CSS, computed by the reverse Kaplan-Meier method.

### Baseline Characteristics and Survival

**Table 1** presents baseline characteristics (full version: **eTable 1**). Of 76,110 patients, 59,821 (78.6%) received no surgery and 16,289 (21.4%) underwent surgical intervention. Non-surgery patients were older (74.6 vs 72.4 years, p<0.001) and more likely Stage IV (26.0% vs 9.7%). Cirrhosis was paradoxically more prevalent in the surgery group (12.9% vs 8.4%, p<0.001).

**Figure 1** displays Kaplan-Meier curves. Median OS by surgery type was: None 4m, Local Destruction 28m, Segmental Resection 30m, Larger Resection 24m, and Transplant 42m (log-rank p<0.0001). In multivariate Cox regression (C-index 0.739; **Figure 2**, **eTable 2**), the strongest protective factors were transplant (HR 0.15 [0.14–0.17]), segmental (HR 0.23 [0.22–0.25]), and larger resection (HR 0.24 [0.22–0.27]). Chemotherapy (HR 0.57 [0.56–0.59]) and radiation (HR 0.52 [0.50–0.53]) were independently protective.

### Primary Finding

**Table 2** presents age-stratified surgery CSS HRs. Across ALL age bands ≥70 years, segmental resection was equivalent or numerically superior to larger resection. At age 70–74: segmental HR 0.27 vs larger HR 0.25. Age 75–79: 0.29 vs 0.33. Age 80+: 0.26 vs 0.31. Only at age 65–69 did larger resection show marginal advantage (0.33 vs 0.35).

**Figure 3** compares HCC and ICC. For HCC: segmental HR 0.26 [0.24–0.28] ≈ larger HR 0.26 [0.22–0.30]. For ICC: segmental HR 0.22 [0.19–0.24] vs larger HR 0.26 [0.22–0.31]. The surgery × cancer type interaction was not significant (p=0.206), indicating that the relative benefit of segmental over larger resection did not differ significantly between HCC and ICC. This non-significant interaction does not confirm equivalence; the study may be underpowered for this test within surgical subgroups.

### Propensity Score Matching

After 1:1 PSM, 6,434 matched pairs were obtained (SMD 0.190→0.030; **eFigure 2**). PSM-adjusted surgery HR was 0.35 (0.33–0.37).

### Machine Learning and External Validation

All survival models demonstrated discrimination in the moderate range: Cox PH (C-index 0.739), RSF (0.736), XGBoost (0.746), and DeepSurv (0.751). The four models were combined into a 5-fold cross-validated ensemble, yielding C-index 0.756 ± 0.003 — a marginal improvement of 0.017 over Cox alone (**eFigure 3**). Model calibration was assessed via time-dependent Brier scores (12/36/60 months): RSF 0.157/0.092/0.059, Cox PH 0.165/0.100/0.066 (**eFigure 9**).

External validation (**Figure 4**, **eFigure 5**): TCGA-LIHC C-indices were 0.595 (95% CI 0.536–0.654, Cox), 0.567 (95% CI 0.508–0.626, RSF), 0.592 (95% CI 0.533–0.651, XGBoost). ICGC C-indices: 0.522 (95% CI 0.461–0.583, Cox), 0.551 (95% CI 0.490–0.612, RSF), 0.547 (95% CI 0.486–0.608, XGBoost). The ΔC≈0.09–0.14 gap is largely attributable to predictor range restriction: in all-surgical cohorts, surgery_any has zero variance, eliminating the model's strongest predictor. The additional 7-point TCGA–ICGC differential (0.595 vs 0.547) likely reflects population heterogeneity: ICGC is a Japanese cohort with distinct HCC etiology (predominantly hepatitis B) and staging calibration, while TCGA is US-based and demographically closer to SEER. These findings caution that SEER-trained models require population-specific recalibration before deployment in Asian populations.

### Sensitivity Analyses

**E-value**: 5.2 for the PSM-adjusted surgery HR of 0.35 (conservative estimate), and 7.5 for the full multivariable model (HR 0.25), with CI-bound values ≥5.8 across all age strata (**eTable 3**). An unmeasured confounder would need to be associated with both surgery receipt and survival by a risk ratio ≥5.2 to nullify the PSM-adjusted effect — exceeding most clinically plausible confounders, including AJCC Stage IV.

**Leave-one-out**: Surgery HR remained 0.27–0.30 excluding any major subgroup (**eTable 5**).

**Model specification**: Surgery HR stable at 0.27→0.29→0.31→0.31→0.25 through sequential adjustment (**eTable 4**).

**Cause-Specific Hazard**: Surgery was protective against cancer-specific death (CSH HR 0.32 [0.31–0.33]) and other-cause death (CSH HR 0.36 [0.34–0.37]), confirming that the benefit extends beyond cancer control. Subdistribution hazard (Fine-Gray) models yielded consistent estimates (data not shown).

**Within-surgery comparison**: For the segmental versus larger resection comparison, the E-value is approximately 1.3, indicating that a weak unmeasured confounder (RR ~1.3) could nullify the observed similarity. This comparison is observational and subject to residual confounding by tumor anatomy and liver function.

**Proportional Hazards**: Schoenfeld residual tests detected violations for 9 of 12 covariates at α=0.05, expected given the large sample size (N=76,110). Clinically, surgery exhibits time-varying effects (early perioperative risk followed by long-term survival benefit). The landmark analysis at 12 and 24 months corroborates the time-averaged HRs. Full Schoenfeld results are provided in eTable 7.

### Temporal, Geographic, and Frailty Analyses

Median OS improved from 3m (2000) to 12m (2019), dipping during COVID-19 (2020–2022: 6m; **eFigure 6**). Higher-income regions had higher surgery rates (24% vs 20%, p<0.001) but all regions showed consistent parenchymal-sparing patterns (**eTable 6**).

Surgery benefit persisted across all GRI strata: Fit HR 0.27, Pre-frail HR 0.32, Frail HR 0.33 (**eFigure 7**). Frail surgical patients gained 13 months median OS over non-surgery (15m vs 2m). RCS identified the inflection point at age 86 (**eFigure 8**).

### Clinical Decision Support

**Table 3** presents the decision matrix. The HBI risk score stratified patients: low-risk (HBI<−5) median OS 22m; intermediate (−5 to 0) 5m; high-risk (>0) 2m (**Figure 5**). The clinical nomogram enables individualized prediction (**Figure 6**).

---

## DISCUSSION

In this population-based analysis of 76,110 elderly hepatobiliary cancer patients, we found that any surgical resection was associated with a large survival benefit that was robust to measured confounding, PSM, and a range of sensitivity analyses. The E-value of 5.2 for the PSM-adjusted estimate indicates that an unmeasured confounder would require RR ≥5.2 to nullify the surgery-versus-none association. This finding confirms and extends prior SEER-based reports [5,6].

The central question — whether segmental/wedge resection provides comparable outcomes to major hepatectomy — is more difficult to answer. The adjusted HRs for segmental (0.23) and larger resection (0.24) were nearly identical, and this pattern was consistent across age bands ≥70, both HCC and ICC, and PSM analyses. However, the E-value for this within-surgery comparison is approximately 1.3, meaning a weak unmeasured confounder (RR ~1.3) could nullify the observed similarity. Such confounders are plausible: tumor proximity to major vasculature, segmental location, underlying liver function (Child-Pugh), and performance status — all absent from SEER — could meaningfully influence both the choice of resection extent and survival. This comparison should be interpreted as an adjusted observational comparison, not evidence of equivalence.

The surgery × cancer type interaction was not significant (p=0.206), indicating that the relative effect of segmental versus larger resection did not differ detectably between HCC and ICC. However, a non-significant interaction does not confirm equivalence, and the study may be underpowered for this test within surgical subgroups. The ICC-specific finding — where the segmental HR (0.22) was numerically lower than larger HR (0.26) — is particularly susceptible to confounding, as ICC surgical planning depends heavily on tumor location and nodal status, data not available in SEER. This finding should not be interpreted as evidence that segmental resection is superior to larger resection for ICC.

The 4-model ensemble achieved C-index 0.756 internally, but this represents a marginal improvement of 0.017 over the Cox model alone. External validation was limited: TCGA-LIHC yielded C-index 0.595, and ICGC-LIRI-JP yielded C-index 0.547 (CI crossing 0.5), indicating that the model does not discriminate survivors from non-survivors in this Japanese cohort. The performance gap is partly attributable to predictor range restriction in all-surgical cohorts (surgery_any has zero variance), but even accounting for this, the models do not provide clinically useful discrimination in external populations. This limitation, combined with the absence of liver function and comorbidity data — variables known to be critical for hepatobiliary cancer prognosis [18,19] — constrains the clinical utility of the prediction models. Chemotherapy (HR 0.57) and radiation (HR 0.52) appeared protective, but this likely reflects selection of healthier patients for treatment rather than true efficacy, and these estimates should be interpreted with caution.

### Strengths and Limitations

Key strengths include the population-based design (N=76,110, the largest such analysis to date), dual external validation cohorts, and multi-layered robustness assessment. Limitations include the absence of liver function variables (Child-Pugh, MELD) — the dominant determinants of surgical candidacy in HCC, where 80–90% of cases arise in cirrhotic livers [18,20]. SEER's cirrhosis field captures only ~9% of the cohort, reflecting >10-fold under-ascertainment, and we do not rely on this variable for confounding control in the primary analysis. The BCLC staging system [18,19] cannot be approximated in SEER. Additional limitations include: the inability to distinguish laparoscopic from open approach; incomplete lymph node data for ICC — a known quality metric [21]; the use of complete-case analysis (eTable 10 reports included vs. excluded comparisons); immortal time bias between diagnosis and intervention (partially mitigated by landmark analyses at 12/24 months); and all-surgical external validation cohorts that eliminate the model's strongest predictor. Temporal trends across the 2004–2022 study period encompassed major changes in HCC management (sorafenib, immunotherapy); diagnosis year was not included as a covariate, which may introduce period effects on the estimated surgery benefit. The Geriatric Risk Index is an exploratory proxy that captures disease severity rather than physiologic frailty. Performance status, comorbidity burden, and resection margin status [22] are absent from all three data sources, representing a universal limitation in population-based hepatobiliary surgical research.

---

## CONCLUSIONS

In this population-based analysis of 76,110 elderly patients with hepatobiliary cancer, any surgical resection was associated with a robust survival benefit (PSM-adjusted HR 0.35; E-value 5.2). Within the operated subgroup, segmental and larger resection showed similar adjusted survival in observational analyses. However, this within-surgery comparison is sensitive to unmeasured confounding (E-value ~1.3) and should not be interpreted as evidence of equivalence. The finding supports parenchymal-sparing approaches when anatomically feasible, but prospective studies with granular clinical data — including liver function, performance status, and tumor anatomy — are needed to inform the choice between segmental and major resection in elderly patients.

---

## REFERENCES

1. Torre LA, Bray F, Siegel RL, Ferlay J, Lortet-Tieulent J, Jemal A. Global cancer statistics, 2012. *CA Cancer J Clin*. 2015;65(2):87–108. doi:10.3322/caac.21262
2. Oweira H, Petrausch U, Helbling D, et al. Early stage hepatocellular carcinoma in the elderly: a SEER database analysis. *J Geriatr Oncol*. 2017;8(4):277–283. doi:10.1016/j.jgo.2017.03.002
3. Ciccarese F, Caturelli E, Felder M, et al. Survival benefit of liver resection for patients with hepatocellular carcinoma across different BCLC stages: a multicentre study. *J Hepatol*. 2014;62(3):617–624. doi:10.1016/j.jhep.2014.10.037
4. Hirokawa F, Hayashi M, Miyamoto Y, et al. Surgical outcomes and clinical characteristics of elderly patients undergoing curative hepatectomy for hepatocellular carcinoma. *J Gastrointest Surg*. 2013;17(11):1929–1937. doi:10.1007/s11605-013-2324-0
5. Zhang QQ, Wu PYS, ALBahde M, et al. Do elderly patients with stage I–II hepatocellular carcinoma benefit from more radical surgeries? A population-based analysis. *Front Oncol*. 2020;10:479. doi:10.3389/fonc.2020.00479
6. Xia L, et al. Effect of surgery on overall survival and cancer-specific survival in patients with primary HCC: a study based on PSM in the SEER cohort. *Medicine*. 2025;104(8):e41521. doi:10.1097/MD.0000000000041521
7. Xu Q, Lu X. Development and validation of an XGBoost model to predict 5-year survival in elderly patients with intrahepatic cholangiocarcinoma after surgery: a SEER-based study. *J Gastrointest Oncol*. 2022;13(6):3290–3299. doi:10.21037/jgo-22-1238
8. Surveillance, Epidemiology, and End Results (SEER) Program. SEER*Stat Database: Incidence – SEER Research Plus Data, 18 Registries, Nov 2024 Sub (2000–2022). National Cancer Institute.
9. Austin PC. Balance diagnostics for comparing the distribution of baseline covariates between treatment groups in propensity-score matched samples. *Stat Med*. 2009;28(25):3083–3107. doi:10.1002/sim.3697
10. Ishwaran H, Kogalur UB, Blackstone EH, Lauer MS. Random survival forests. *Ann Appl Stat*. 2008;2(3):841–860. doi:10.1214/08-AOAS169
11. Chen T, Guestrin C. XGBoost: a scalable tree boosting system. *KDD*. 2016. doi:10.1145/2939672.2939785
12. Uno H, Cai T, Pencina MJ, D'Agostino RB, Wei LJ. On the C-statistics for evaluating overall adequacy of risk prediction procedures with censored survival data. *Stat Med*. 2011;30(10):1105–1117. doi:10.1002/sim.4154
13. Steyerberg EW, Vickers AJ, Cook NR, et al. Assessing the performance of prediction models: a framework for traditional and novel measures. *Epidemiology*. 2010;21(1):128–138. doi:10.1097/EDE.0b013e3181c30fb2
14. Cancer Genome Atlas Research Network. Comprehensive and integrative genomic characterization of hepatocellular carcinoma. *Cell*. 2017;169(7):1327–1341. doi:10.1016/j.cell.2017.05.046
15. Totoki Y, Tatsuno K, Covington KR, et al. Trans-ancestry mutational landscape of hepatocellular carcinoma genomes. *Nat Genet*. 2014;46:1267–1273. doi:10.1038/ng.3126
16. VanderWeele TJ, Ding P. Sensitivity analysis in observational research: introducing the E-value. *Ann Intern Med*. 2017;167(4):268–274. doi:10.7326/M16-2607
17. Fine JP, Gray RJ. A proportional hazards model for the subdistribution of a competing risk. *J Am Stat Assoc*. 1999;94(446):496–509. doi:10.1080/01621459.1999.10474144
18. European Association for the Study of the Liver. EASL Clinical Practice Guidelines: management of hepatocellular carcinoma. *J Hepatol*. 2018;69(1):182–236. doi:10.1016/j.jhep.2018.03.019
19. Heimbach JK, Kulik LM, Finn RS, et al. AASLD guidelines for the treatment of hepatocellular carcinoma. *Hepatology*. 2018;67(1):358–380. doi:10.1002/hep.29086
20. Marrero JA, Kulik LM, Sirlin CB, et al. Diagnosis, staging, and management of hepatocellular carcinoma: 2018 practice guidance by the American Association for the Study of Liver Diseases. *Hepatology*. 2018;68(2):723–750. doi:10.1002/hep.29913
21. Mavros MN, Mayo SC, Hyder O, et al. A systematic review of the outcomes of hepatic resection for intrahepatic cholangiocarcinoma. *Ann Surg Oncol*. 2014;21(8):2672–2680. doi:10.1245/s10434-014-3608-8
22. Cauchy F, Fuks D, Nomi T, et al. Incidence, risk factors and consequences of positive resection margin after hepatectomy for colorectal liver metastases. *Eur J Surg Oncol*. 2016;42(3):357–364. doi:10.1016/j.ejso.2015.12.008

---

## SUPPLEMENTARY MATERIAL

**eFigure 1–9**: CONSORT Diagram, PSM Love Plot, 5-Fold CV, SHAP Analysis, External Validation Detail, Temporal Trends + COVID, Geriatric Risk Index, Age-Surgery RCS + Competing Risk CIF, Model Calibration

**eTable 1–10**: Full Baseline Characteristics, Complete Multivariate Cox, Stratified E-values, Model Specification Robustness, Leave-One-Out Analysis, Geographic Practice Variation, Model Hyperparameters, Schoenfeld + Landmark, Clinical Decision Matrix, Included vs. Excluded Patient Comparison

---

*Correspondence to:* Zhuha Zhou, MD, Department of Gastroenterology Surgery, The First Affiliated Hospital of Wenzhou Medical University, Zhejiang, China. E-mail: zhouzhuha@wmu.edu.cn

**Author Contributions**: QX and YB contributed equally (co-first authors). ZZ had full access to all data and takes responsibility for integrity. Study concept and design: QX, ZZ. Data acquisition and analysis: QX, YB, YC. Statistical analysis: ZZ. Machine learning: ZZ, YB. Manuscript drafting: ZZ. Critical revision: All authors. Study supervision: ZZ.

**Funding**: This research did not receive any specific grant from funding agencies in the public, commercial, or not-for-profit sectors.

**Conflict of Interest**: The authors have declared that no competing interest exists.

**Ethics Statement**: This study used publicly available, de-identified data from SEER, TCGA, and ICGC and was exempt from institutional review board approval.

**Data Availability**: SEER data: seer.cancer.gov (Research Plus Data, November 2024 submission). TCGA-LIHC: portal.gdc.cancer.gov. ICGC-LIRI-JP: dcc.icgc.org. Analysis code: github.com/woodhaha/hepatobiliary-cancer-seer-analysis.
