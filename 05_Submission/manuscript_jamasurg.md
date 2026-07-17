# Comparison of Segmental Resection vs Major Hepatectomy in Elderly Hepatobiliary Cancer

**Authors**: Qigang Xu, MD$^{1,2,\dagger}$; Yongyu Bai, MD$^{1,\dagger}$; Yiqi Cai, MD$^{1}$; Zhuha Zhou, MD$^{1,*}$

$^{1}$ Department of Gastroenterology Surgery, The First Affiliated Hospital of Wenzhou Medical University, Zhejiang, China  
$^{2}$ Department of Hepatobiliary and Pancreatic Surgery, The First Affiliated Hospital of Wenzhou Medical University, Zhejiang, China  
$\dagger$ Equal contribution  
$^{*}$ **Corresponding author**: Zhuha Zhou, MD, Department of Gastroenterology Surgery, The First Affiliated Hospital of Wenzhou Medical University, Wenzhou, Zhejiang 325000, China. Email: zhouzhuha@wmu.edu.cn. Phone: +86-13806678098.

**Manuscript word count**: 2,847  
**Key Points word count**: 96  
**Abstract word count**: 348

---

## Key Points

**Question**: For elderly patients with hepatobiliary cancer, does segmental resection provide survival comparable to major hepatectomy?

**Findings**: In this cohort study of 76,110 elderly patients, major hepatectomy did not confer better cancer-specific survival than segmental/wedge resection across age bands ≥70 and both HCC and ICC (HR 0.24 vs 0.23), though findings are sensitive to unmeasured confounding (E-value ≈ 1.3). A 4-model ensemble achieved C-index 0.756.

**Meaning**: When anatomically feasible, segmental resection may be preferred to reduce perioperative risk in elderly patients, but clinical decisions must incorporate tumor anatomy, liver function, and comorbidity. Major hepatectomy was not associated with better survival.

---

## Abstract

**Importance**: The choice between parenchymal-sparing segmental resection and major hepatectomy in elderly patients with hepatobiliary cancer remains unresolved, yet directly affects surgical decision-making in a rapidly aging population.

**Objective**: To compare cancer-specific survival between segmental/wedge resection and major hepatectomy in patients aged 65 years or older with hepatocellular carcinoma (HCC) or intrahepatic cholangiocarcinoma (ICC), and to quantify the robustness of findings to unmeasured confounding.

**Design, Setting, and Participants**: Population-based retrospective cohort study using the SEER 18-registry database (2004–2022) with external validation in TCGA-LIHC and ICGC-LIRI-JP. Of 171,286 records, 76,110 patients aged ≥65 years with histologically confirmed hepatobiliary malignancy met inclusion criteria. Data were analyzed April to June 2026.

**Exposures**: Five surgical categories — no surgery, local destruction, segmental/wedge resection, lobectomy/extended lobectomy, and liver transplantation — classified using SEER site-specific codes.

**Main Outcomes and Measures**: Cancer-specific survival. Multivariable Cox proportional hazards regression with propensity score matching (1:1, caliper 0.05) and E-value analysis for unmeasured confounding. Secondary analyses included a 4-model machine learning ensemble (Cox, Random Survival Forest, XGBoost, DeepSurv) with 5-fold cross-validation, competing risk analysis, and landmark conditional survival.

**Results**: Among 76,110 patients (median age 74 years [IQR 68–80]; 49,774 male [65.4%]), 59,821 (78.6%) had no surgery and 16,289 (21.4%) had surgery. Surgical resection was associated with improved cancer-specific survival vs no surgery (PSM-adjusted HR 0.35; 95% CI 0.33–0.37; E-value 5.2). Among surgical patients, segmental resection (HR 0.23; 95% CI 0.22–0.25) and major hepatectomy (HR 0.24; 95% CI 0.22–0.27) were similar across age bands ≥70 years. The surgery × cancer type interaction was not significant (P = 0.21). This comparison is susceptible to unmeasured confounding (E-value ≈ 1.3). External validation showed moderate discrimination: TCGA-LIHC C-index 0.595 (95% CI 0.536–0.654); ICGC-LIRI-JP C-index 0.547 (95% CI 0.486–0.608).

**Conclusions and Relevance**: Surgical resection is associated with a robust survival benefit in elderly hepatobiliary cancer patients. Although adjusted outcomes for segmental and major resection appear similar, this finding should be interpreted cautiously given susceptibility to unmeasured confounding. When oncologically appropriate, parenchymal-sparing resection may be considered in elderly patients, but prospective studies with granular clinical data are needed to guide the extent of resection.

---

## Introduction

Hepatobiliary cancers (hepatocellular carcinoma [HCC] and intrahepatic cholangiocarcinoma [ICC]) are among the most lethal malignancies worldwide, and their incidence in older adults is rising with population aging.[1,2] For elderly patients with resectable disease, surgeons choose between parenchymal-sparing approaches (segmental or wedge resection) and major hepatectomy (lobectomy or extended lobectomy). Major resection may offer superior oncologic clearance but carries higher perioperative risk, particularly in elderly patients with comorbidities and reduced physiologic reserve.[3,4] A recent systematic review confirmed that surgical resection in elderly patients with HCC is safe and achieves oncologic outcomes comparable to younger patients,[21] yet the optimal extent of resection remains undefined. Current guidelines recommend resection based on tumor stage and liver function[5,6] but provide limited evidence to guide the choice of resection extent in older patients.

Zhang et al. reported from SEER data that for 2,371 elderly patients with stage I–II HCC, segmental or wedge resection provided survival comparable to larger resection when transplantation was unavailable.[7] However, that analysis was limited to early-stage HCC, used data through 2013, and did not include ICC. Xia et al. compared surgical types in 10,174 patients with HCC but combined wedge and segmental resection as a single category without formal causal sensitivity analysis.[8] The present study addresses these gaps through three contributions: (1) simultaneous evaluation of both HCC and ICC in the largest elderly surgical cohort to date (N = 76,110); (2) formal E-value analysis to quantify robustness to unmeasured confounding across all comparisons; and (3) external validation of a multi-model survival ensemble in independent genomic cohorts.

To this end, we conducted a population-based analysis with three objectives: (1) to quantify the survival benefit of surgical resection in elderly hepatobiliary cancer patients using propensity score matching and E-value analysis; (2) to compare adjusted outcomes by surgical extent within the operated subgroup across HCC and ICC; and (3) to assess robustness through multiple sensitivity analyses including external validation in independent cohorts.

---

## Methods

### Data Source and Study Population

The SEER Research Plus Data (November 2024 submission, 18 registries, covering approximately 28% of the US population) was queried for patients with hepatobiliary cancer (ICD-O-3 site codes: C22.0 liver, C22.1 intrahepatic bile duct) diagnosed from January 1, 2004 to December 31, 2022.[9] Inclusion criteria were: (1) age ≥65 years at diagnosis; (2) histologically confirmed hepatobiliary malignancy (HCC: ICD-O-3 8170–8175, 8180; ICC: 8160–8162); (3) known AJCC stage and survival data; (4) complete data on key covariates. Of 171,286 records, 76,110 (44.4%) met criteria; exclusions included age <65 (n = 89,172), missing survival (n = 5,504), and incomplete covariates (n = 500, of which the most frequent were missing tumor size [n = 312] and unknown grade [n = 188]). Complete-case analysis was used given the low overall missing rate (0.7% of age-eligible records). External validation was performed on TCGA-LIHC (n = 269)[10] and ICGC-LIRI-JP (n = 260).[11] The study was exempt from institutional review board approval as all data sources are publicly available and de-identified.

### Variable Definitions

The primary exposure was surgery categorized into five groups based on SEER site-specific codes: (1) None (no cancer-directed surgery); (2) Local Destruction (codes 10–27; radiofrequency ablation, percutaneous ethanol injection, cryoablation); (3) Segmental Resection (codes 30–38; wedge or segmental resection); (4) Major Resection (codes 50–59; lobectomy, extended lobectomy); (5) Liver Transplantation (codes 60–66). A unified AJCC staging algorithm harmonized the 3rd, 6th, 7th editions and EOD 2018 staging systems to stage groups I–IV, consistent with prior SEER hepatobiliary studies.[7,8] The primary outcome was cancer-specific survival (CSS). Covariates included age, sex, marital status, race, AJCC stage, tumor grade, cancer type, chemotherapy, radiation, cirrhosis, and tumor size.

### Statistical Analysis

This study was designed as an emulation of a target trial comparing segmental resection versus major hepatectomy for resectable hepatobiliary cancer in patients aged 65 years or older, following the TARGET framework proposed by Cashin et al.[22] The protocol specified: (1) eligibility criteria — age ≥65 years, histologically confirmed HCC or ICC, known AJCC stage and survival data, complete covariate data; (2) treatment strategies — segmental/wedge resection (SEER codes 30–38) versus major hepatectomy (codes 50–59); (3) assignment procedure — treatment assignment was based on clinical practice (observational); (4) follow-up period — from diagnosis to death, censoring, or end of study (December 2022); (5) causal contrast — per-protocol effect (as-treated analysis); (6) statistical analysis plan — multivariable Cox regression, propensity score matching, and E-value analysis as prespecified in the supplementary protocol (eTable 9, TARGET Protocol Specification).

Baseline characteristics were compared using standardized mean differences (SMD). Cancer-specific survival was estimated using the Kaplan-Meier method with log-rank testing. Multivariable Cox proportional hazards regression was performed for CSS, with stratified analyses by age band (65–69, 70–74, 75–79, 80+), AJCC stage (I–IV), and cancer type (HCC, ICC). Formal interaction testing between surgery and cancer type, age, and stage was conducted. The proportional hazards assumption was evaluated using Schoenfeld residuals.

With N = 76,110 and approximately 82.7% mortality, the study had >99% power at α = 0.05 to detect a hazard ratio of 1.1 or smaller for the primary surgery-versus-none comparison. For the segmental-versus-major comparison (n ≈ 3,185), statistical power was 80% to detect a HR of 1.25 or greater.

Propensity score matching (1:1, nearest-neighbor, caliper 0.05, logistic regression) was performed with covariates including age, sex, marital status, stage, grade, cancer type, chemotherapy, radiation, and tumor size. Balance was assessed using SMD (<0.1 considered adequate). A 4-model survival ensemble (L2-regularized Cox PH, Random Survival Forest [100 trees],[12] XGBoost [survival:cox, max_depth=4, learning_rate=0.05],[13] and DeepSurv [3-layer MLP: 17→64→32→16→1, ReLU, batch normalization, dropout=0.3])[14] was trained using temporal split (2004–2017 training, 2018–2022 testing) with 5-fold cross-validation. Training and test set C-indices are reported separately to assess overfitting (eTable 8). Model performance was assessed using Harrell's C-index with Uno's correction[15] and time-dependent Brier scores at 12, 36, and 60 months.

### Sensitivity Analyses

E-value analysis quantified the minimum strength of association between an unmeasured confounder and both treatment and outcome needed to nullify the observed association,[16,23] following established guidance for applying E-values in clinical research.[23] Additional sensitivity analyses included: cause-specific hazard and Fine-Gray competing risk models,[17] leave-one-out subgroup exclusion, sequential model specification testing, restricted cubic splines for non-linear age effects, Geriatric Risk Index (GRI) stratification, and landmark conditional survival analysis at 12 and 24 months. The Hepatobiliary Cancer Index (HBI), a 7-variable risk score, was derived from the ensemble model and validated internally.

### External Validation

SEER-trained models were applied to TCGA-LIHC (all-surgical HCC)[10] and ICGC-LIRI-JP (all-surgical HCC, Japanese cohort)[11] without retraining. Model discrimination was assessed using Harrell's C-index. The study follows TRIPOD Type 3 reporting guidelines. All tests were two-sided with α = 0.05. Analyses were performed using Python 3.13 (scikit-survival v0.22, lifelines v0.28, xgboost v2.1, PyTorch 2.5) and R 4.4.

---

## Results

### Study Population

The analytic cohort comprised 76,110 elderly patients: 57,380 HCC (75.4%), 11,749 ICC (15.4%), and 6,981 other hepatobiliary malignancies (9.2%; including combined HCC-ICC, hepatoblastoma, sarcoma NOS). Median age was 74 years (IQR 68–80); 49,774 (65.4%) were male. Median follow-up was 49 months (IQR 16–96) for overall survival and 74 months (IQR 24–120) for CSS. The cohort selection process is summarized in eFigure 1.

Of 76,110 patients, 59,821 (78.6%) received no surgery and 16,289 (21.4%) underwent surgical intervention (Table 1). Non-surgery patients were older (mean 74.6 vs 72.4 years; SMD 0.25) and more likely to have Stage IV disease (26.0% vs 9.7%; SMD 0.42). Cirrhosis was more prevalent in the surgery group (12.9% vs 8.4%; SMD 0.15), reflecting selection of patients with preserved liver function for intervention.

### Surgical Resection and Survival

Surgical resection was associated with improved cancer-specific survival compared to no surgery (PSM-adjusted HR 0.35; 95% CI 0.33–0.37; E-value 5.2). This benefit was consistent across age, stage, and cancer type strata. After 1:1 propensity score matching, 6,434 matched pairs were obtained with adequate balance (SMD 0.190 → 0.030). An unmeasured confounder would need a risk ratio of at least 5.2 to explain this association — exceeding most clinically plausible confounders including AJCC Stage IV.

Within the surgical subgroup, segmental resection (HR 0.23; 95% CI 0.22–0.25) and major hepatectomy (HR 0.24; 95% CI 0.22–0.27) showed remarkably similar adjusted survival in multivariable analysis (Table 2). This pattern persisted across all age bands ≥70 years and both HCC and ICC (Figure 1, Figure 2). At age 70–74: segmental HR 0.27 vs major HR 0.25; age 75–79: 0.29 vs 0.33; age 80+: 0.26 vs 0.31. For HCC, segmental HR 0.26 (95% CI 0.24–0.28) and major HR 0.26 (95% CI 0.22–0.30); for ICC, segmental HR 0.22 (95% CI 0.19–0.24) and major HR 0.26 (95% CI 0.22–0.31). The surgery × cancer type interaction was not significant (P = 0.21), and the surgery × age band interaction was likewise not significant (P = 0.38). The E-value for the segmental-vs-major comparison was approximately 1.3, indicating that a weak unmeasured confounder (risk ratio ≈ 1.3) could plausibly nullify the observed similarity. Clinically relevant confounders including Child-Pugh grade, performance status, tumor proximity to vasculature, and resection margin status — all absent from SEER — may influence both the choice of resection extent and survival outcomes.

### Machine Learning and External Validation

The 4-model survival ensemble achieved an internal C-index of 0.756 (95% CI 0.753–0.759), with individual models ranging from 0.739 (Cox) to 0.751 (DeepSurv). Training-test C-index convergence was similar across models (ΔC between training and test folds < 0.02 for all models; eTable 8), indicating no substantial overfitting despite the DeepSurv architecture's parameter count. Model calibration at 12, 36, and 60 months was satisfactory: RSF Brier scores of 0.157, 0.092, and 0.059 respectively. External validation in all-surgical HCC cohorts yielded moderate discrimination: TCGA-LIHC C-index 0.595 (95% CI 0.536–0.654) and ICGC-LIRI-JP C-index 0.547 (95% CI 0.486–0.608). This attenuation relative to internal validation is expected given predictor range restriction — surgery type has limited variance in all-surgical cohorts — and population heterogeneity between the predominantly Western SEER cohort and the Japanese ICGC cohort with distinct HCC etiology.[11] The moderate external performance highlights the need for population-specific recalibration before clinical application.

### Sensitivity Analyses

The E-value bound of 5.2 remained stable across all age strata (range 4.8–6.1). Cause-specific hazard models confirmed the protective effect of surgery against both cancer-specific death (HR 0.32; 95% CI 0.31–0.33) and other-cause death (HR 0.36; 95% CI 0.34–0.37). Sequential model specification demonstrated stable surgery HR estimates (0.27 → 0.29 → 0.31 → 0.31 → 0.25) through progressive covariate adjustment (eTable 3). In cancer-type stratified analyses, the pattern of similar adjusted survival between segmental and major resection was observed in both HCC (segmental HR 0.26; 95% CI 0.24–0.28; major HR 0.26; 95% CI 0.22–0.30) and ICC (segmental HR 0.22; 95% CI 0.19–0.24; major HR 0.26; 95% CI 0.22–0.31), consistent with the non-significant interaction (P = 0.21). The surgery benefit persisted across all Geriatric Risk Index strata (fit HR 0.27, pre-frail HR 0.32, frail HR 0.33; eFigure 5), and the HBI risk score identified distinct prognostic groups in internal validation (Figure 3). Exploratory analysis of patients aged ≥85 years (n = 276 surgical patients) showed similar patterns, though limited by sample size (eTable 10).

---

## Discussion

In this population-based analysis of 76,110 elderly patients with hepatobiliary cancer, surgical resection was associated with a substantial and robust survival benefit (PSM-adjusted HR 0.35; E-value 5.2), confirming prior reports[7,8] and extending them with formal causal sensitivity analysis. The more clinically relevant question is whether major hepatectomy offers better cancer-specific survival than segmental resection in elderly patients. In adjusted analyses, major hepatectomy was not associated with improved survival compared with segmental resection (HR 0.24 vs 0.23), across all age bands ≥70 years and both HCC and ICC. However, the E-value of approximately 1.3 signals that a weak unmeasured confounder could fully account for this similarity, and the finding should not be interpreted as evidence of equivalence. Clinically, when oncologically appropriate, segmental resection may be preferred to limit perioperative risk. Prior NSQIP-based studies have demonstrated that major hepatectomy in elderly patients carries substantially higher morbidity — severe complication-associated mortality of 20.1% in elderly versus 10.8% in non-elderly patients[28] — and existing risk calculators may underestimate these risks.[29] These perioperative risks, combined with the absence of a survival advantage for major resection observed in the present study, support consideration of parenchymal-sparing approaches when anatomically feasible. This is consistent with accumulating evidence that major resection does not improve oncologic outcomes in elderly patients when parenchymal-sparing alternatives are anatomically feasible.[18,19] Jung et al. recently reported that non-anatomic resection for HCC did not compromise recurrence or survival compared to anatomic resection in a propensity-matched analysis,[24] and Weigle et al. found that octogenarians undergoing resection for ICC received fewer major hepatectomies with comparable 1-year survival to younger patients.[25]

The 4-model ensemble incorporating clinical and demographic variables achieved moderate internal discrimination (C-index 0.756), with limited added value from complex machine learning methods beyond Cox regression (ΔC = 0.017). External validation performance was moderate, with a notable gap between TCGA (C-index 0.595) and ICGC (C-index 0.547). This attenuation is expected for several reasons. First, both external cohorts are composed exclusively of surgical patients, substantially limiting variance in the most influential predictor — surgery type. Second, the predictor sets differ between SEER (clinical and demographic variables) and the genomic cohorts, where many SEER-available variables are not recorded. Third, population heterogeneity — particularly between the predominantly Western SEER cohort and the Japanese ICGC cohort with distinct HCC etiology (viral hepatitis predominance)[11] — further attenuates transportability. These factors are well-recognized challenges in cross-cohort model validation[15] and do not undermine the internal validity of the primary survival comparisons. The HBI risk score, derived from the ensemble model, stratified patients into distinct prognostic groups in internal validation and warrants further external validation as an exploratory tool to inform surgical decision-making.

### Limitations

This study has several important limitations. First, liver function data (Child-Pugh class, MELD score) and performance status are not available in SEER. Recent evidence indicates that frailty is a strong and independent predictor of postoperative outcomes after hepatic resection for HCC, independent of chronologic age,[26] and frailty assessment is increasingly recognized as essential in hepatobiliary surgical decision-making.[27] The absence of these variables from SEER likely contributes to residual confounding. These factors are critical determinants of surgical candidacy in HCC, where most cases arise in cirrhotic livers.[5,6] Their absence limits confounding control for the surgery-vs-none comparison and precludes definitive assessment of the segmental-vs-major comparison. The E-value of 1.3 quantifies this vulnerability: any confounder associated with a risk ratio of approximately 1.3 with both treatment selection and survival — a threshold easily crossed by performance status or tumor anatomy — could nullify the apparent similarity. Second, BCLC staging — the standard HCC treatment allocation framework — cannot be approximated from SEER variables, limiting the granularity of treatment context. Third, margin status, laparoscopic vs open approach, and detailed tumor anatomy are unavailable. Fourth, SEER does not capture perioperative outcomes such as 90-day mortality, major complications, or length of stay, which are critical for weighing the trade-off between segmental and major resection in elderly patients; prior NSQIP-based studies have demonstrated that major hepatectomy carries substantially higher morbidity in older adults (severe complication-associated mortality 20.1% in elderly vs 10.8% in non-elderly),[28] and existing risk calculators may underestimate these risks.[29] Future studies using NSQIP or institutional data should address this gap. Fifth, both external validation cohorts are all-surgical HCC samples, limiting generalizability to ICC or unresected patients. Dedicated ICC genomic cohorts suitable for survival model validation remain limited; future studies incorporating ICC-specific external validation (e.g., from large multicenter cohorts or upcoming ICC genomic initiatives) are needed. Sixth, temporal changes in systemic therapy (including immune checkpoint inhibitors introduced after 2017) may introduce period effects not fully captured by the analysis. Seventh, E-value analysis does not address measurement error or model misspecification; unmeasured confounding remains possible even with high E-values.[20]

### Conclusions

Surgical resection was associated with a large and robust survival benefit in elderly patients with hepatobiliary cancer. In adjusted analyses, major hepatectomy was not associated with improved cancer-specific survival compared with segmental resection, though this finding is sensitive to unmeasured confounding (E-value ≈ 1.3). When oncologically appropriate, parenchymal-sparing resection may be preferred to limit perioperative risk, but clinical decisions must consider tumor anatomy, liver function, and comorbidity. Prospective studies with comprehensive clinical data are needed to guide resection extent in this growing population.

---

## References

1. Torre LA, Bray F, Siegel RL, Ferlay J, Lortet-Tieulent J, Jemal A. Global cancer statistics, 2012. *CA Cancer J Clin*. 2015;65(2):87–108. doi:10.3322/caac.21262
2. Oweira H, Petrausch U, Helbling D, et al. Early stage hepatocellular carcinoma in the elderly: a SEER database analysis. *J Geriatr Oncol*. 2017;8(4):277–283. doi:10.1016/j.jgo.2017.03.002
3. Hirokawa F, Hayashi M, Miyamoto Y, et al. Surgical outcomes and clinical characteristics of elderly patients undergoing curative hepatectomy for hepatocellular carcinoma. *J Gastrointest Surg*. 2013;17(11):1929–1937. doi:10.1007/s11605-013-2324-0
4. Cauchy F, Fuks D, Nomi T, et al. Incidence, risk factors and consequences of positive resection margin after hepatectomy for colorectal liver metastases. *Eur J Surg Oncol*. 2016;42(3):357–364. doi:10.1016/j.ejso.2015.12.008
5. European Association for the Study of the Liver. EASL Clinical Practice Guidelines: management of hepatocellular carcinoma. *J Hepatol*. 2018;69(1):182–236. doi:10.1016/j.jhep.2018.03.019
6. Marrero JA, Kulik LM, Sirlin CB, et al. Diagnosis, staging, and management of hepatocellular carcinoma: 2018 practice guidance by the American Association for the Study of Liver Diseases. *Hepatology*. 2018;68(2):723–750. doi:10.1002/hep.29913
7. Zhang QQ, Wu PYS, ALBahde M, et al. Do elderly patients with stage I–II hepatocellular carcinoma benefit from more radical surgeries? A population-based analysis. *Front Oncol*. 2020;10:479. doi:10.3389/fonc.2020.00479
8. Xia L, et al. Effect of surgery on overall survival and cancer-specific survival in patients with primary HCC: a study based on PSM in the SEER cohort. *Medicine*. 2025;104(8):e41521. doi:10.1097/MD.0000000000041521
9. Surveillance, Epidemiology, and End Results (SEER) Program. SEER*Stat Database: Incidence – SEER Research Plus Data, 18 Registries, Nov 2024 Sub (2000–2022). National Cancer Institute.
10. Cancer Genome Atlas Research Network. Comprehensive and integrative genomic characterization of hepatocellular carcinoma. *Cell*. 2017;169(7):1327–1341. doi:10.1016/j.cell.2017.05.046
11. Totoki Y, Tatsuno K, Covington KR, et al. Trans-ancestry mutational landscape of hepatocellular carcinoma genomes. *Nat Genet*. 2014;46:1267–1273. doi:10.1038/ng.3126
12. Ishwaran H, Kogalur UB, Blackstone EH, Lauer MS. Random survival forests. *Ann Appl Stat*. 2008;2(3):841–860. doi:10.1214/08-AOAS169
13. Chen T, Guestrin C. XGBoost: a scalable tree boosting system. *KDD*. 2016. doi:10.1145/2939672.2939785
14. Katzman JL, Shaham U, Cloninger A, Bates J, Jiang T, Kluger Y. DeepSurv: personalized treatment recommender system using a Cox proportional hazards deep neural network. *BMC Med Res Methodol*. 2018;18(1):24. doi:10.1186/s12874-018-0482-1
15. Uno H, Cai T, Pencina MJ, D'Agostino RB, Wei LJ. On the C-statistics for evaluating overall adequacy of risk prediction procedures with censored survival data. *Stat Med*. 2011;30(10):1105–1117. doi:10.1002/sim.4154
16. VanderWeele TJ, Ding P. Sensitivity analysis in observational research: introducing the E-value. *Ann Intern Med*. 2017;167(4):268–274. doi:10.7326/M16-2607
17. Fine JP, Gray RJ. A proportional hazards model for the subdistribution of a competing risk. *J Am Stat Assoc*. 1999;94(446):496–509. doi:10.1080/01621459.1999.10474144
18. Mavros MN, Mayo SC, Hyder O, et al. A systematic review of the outcomes of hepatic resection for intrahepatic cholangiocarcinoma. *Ann Surg Oncol*. 2014;21(8):2672–2680. doi:10.1245/s10434-014-3608-8
19. Egger ME, Ohlson EC, Scoggins CR, et al. Assessment of chemotherapy response in colorectal liver metastases in patients undergoing hepatic resection and the correlation with pathologic residual tumor size. *J Am Coll Surg*. 2018;226(5):845–852. doi:10.1016/j.jamcollsurg.2017.12.040
20. VanderWeele TJ, Mathur MB, Chen Y. Outcome-wide longitudinal designs for causal inference: a new template for empirical studies. *Stat Sci*. 2020;35(3):437–466. doi:10.1214/19-STS728
21. Lee JS, Park DA, Yoo JJ, et al. Efficacy and safety of surgical resection in elderly patients with hepatocellular carcinoma: a systematic review and meta-analysis. *Gut Liver*. 2024;18(4):695–708. doi:10.5009/gnl230485
22. Cashin AG, et al. Transparent reporting of observational studies emulating a target trial (TARGET). *JAMA*. 2025;333(20):1845–1854. doi:10.1001/jama.2025.0483
23. Haneuse S, VanderWeele TJ, Arterburn D. Using the E-value to assess the potential effect of unmeasured confounding in observational studies. *JAMA*. 2019;321(6):602–603. doi:10.1001/jama.2018.21554
24. Jung S, et al. Parenchymal-sparing non-anatomic resection vs. classic anatomic resection for hepatocellular carcinoma: a propensity score matched analysis. *Ann Hepatobiliary Pancreat Surg*. 2025;29(2):121–130.
25. Weigle CA, et al. Resection of intrahepatic cholangiocarcinoma in octogenarians. *HPB*. 2024;26(8):1037–1045.
26. Diao YK, et al. Association of preoperative frailty with short- and long-term outcomes after hepatic resection for hepatocellular carcinoma. *BJS Open*. 2025;9(1):zrae171.
27. Funamizu N, et al. Frailty in hepatobiliary and pancreatic surgery: a narrative review. *Gland Surg*. 2026.
28. Tzeng CWD, et al. Predictors of morbidity and mortality after hepatectomy in elderly patients: analysis of 7,621 NSQIP patients. *HPB*. 2014;16(5):459–466.
29. Sahara K, et al. Evaluation of the ACS NSQIP surgical risk calculator in elderly patients undergoing hepatectomy for hepatocellular carcinoma. *J Gastrointest Surg*. 2020;24(3):551–558.

---

## Figure Legends

**Figure 1** — Cancer-Specific Survival by Surgery Type Among Elderly Patients With Hepatobiliary Cancer (N = 76,110). Kaplan-Meier curves comparing cancer-specific survival across five surgical categories: no surgery, local destruction, segmental/wedge resection, major hepatectomy (lobectomy/extended lobectomy), and liver transplantation. Shaded regions indicate 95% CIs. Log-rank P < 0.001. Tick marks indicate censored observations. Number at risk shown at 0, 12, 24, 36, 48, and 60 months.

**Figure 2** — Forest Plot of Multivariable Cox Regression Results for Cancer-Specific Survival. Adjusted hazard ratios for surgical resection (segmental vs major) across subgroups defined by age, sex, stage, cancer type, and frailty. Models adjusted for age, sex, marital status, race, AJCC stage, tumor grade, cancer type, chemotherapy, radiation, cirrhosis, and tumor size. Squares represent point estimates; horizontal lines indicate 95% CIs.

**Figure 3** — Hepatobiliary Cancer Index (HBI) Risk Score and Clinical Decision Support. A) Distribution of HBI risk scores in the study population. B) Cancer-specific survival stratified by HBI risk group (low, intermediate, high). C) Decision matrix integrating HBI score with age and stage. D) Nomogram for individualized survival prediction. The HBI score was derived from the 4-model ensemble using 7 clinical variables.

---

## Tables

### Table 1. Baseline Characteristics of Elderly Patients With Hepatobiliary Cancer by Surgery Status

| Characteristic | No Surgery (n = 59,821) | Surgery (n = 16,289) | SMD |
|---|---|---|---|
| Age, mean (SD), y | 74.6 (7.8) | 72.4 (7.2) | 0.25 |
| Male, No. (%) | 39,284 (65.7) | 10,490 (64.4) | 0.03 |
| AJCC Stage, No. (%) | | | 0.42 |
| — I | 21,865 (36.6) | 8,941 (54.9) | |
| — II | 12,490 (20.9) | 3,429 (21.0) | |
| — III | 9,915 (16.6) | 2,326 (14.3) | |
| — IV | 15,551 (26.0) | 1,593 (9.8) | |
| HCC, No. (%) | 43,921 (73.4) | 13,459 (82.6) | 0.22 |
| ICC, No. (%) | 9,553 (16.0) | 2,196 (13.5) | 0.07 |
| Cirrhosis, No. (%) | 5,014 (8.4) | 2,096 (12.9) | 0.15 |
| Chemotherapy, No. (%) | 18,543 (31.0) | 3,586 (22.0) | 0.20 |
| Radiation, No. (%) | 7,866 (13.1) | 931 (5.7) | 0.26 |
| Tumor size, mean (SD), cm | 5.7 (3.8) | 4.2 (3.1) | 0.33 |

SMD indicates standardized mean difference; HCC, hepatocellular carcinoma; ICC, intrahepatic cholangiocarcinoma. Complete baseline characteristics including age bands, race, marital status, grade, and geographic region are presented in eTable 1.

### Table 2. Multivariable Cox Regression for Cancer-Specific Survival

| Variable | HR (95% CI) | P value |
|---|---|---|
| **Surgery (vs None)** | | |
| — Local Destruction | 0.53 (0.51–0.55) | <0.001 |
| — Segmental Resection | 0.23 (0.22–0.25) | <0.001 |
| — Major Hepatectomy | 0.24 (0.22–0.27) | <0.001 |
| — Liver Transplantation | 0.15 (0.14–0.17) | <0.001 |
| **Age** (per year) | 1.03 (1.03–1.03) | <0.001 |
| **Sex** (male vs female) | 1.05 (1.03–1.07) | <0.001 |
| **Cancer Type** (ICC vs HCC) | 1.12 (1.09–1.15) | <0.001 |
| **AJCC Stage** (vs I) | | |
| — II | 1.27 (1.24–1.30) | <0.001 |
| — III | 1.69 (1.65–1.73) | <0.001 |
| — IV | 2.36 (2.30–2.42) | <0.001 |
| **Chemotherapy** | 0.57 (0.56–0.59) | <0.001 |
| **Cirrhosis** | 0.88 (0.85–0.91) | <0.001 |

HR indicates hazard ratio; CI, confidence interval. Model adjusted for age, sex, marital status, race, AJCC stage, grade, cancer type, surgery, chemotherapy, radiation, cirrhosis, and tumor size. Complete regression results including all covariates are provided in eTable 2.

---

## Supplementary Material (Online-Only)

**eFigures**: eFig1 CONSORT Diagram · eFig2 SHAP Feature Importance · eFig3 Model Comparison · eFig4 Temporal Trends · eFig5 Landmark + Frailty · eFig6 Calibration · eFig7 Subgroup Forest · eFig8 Overview

**eTables**: eTable1 Full Baseline · eTable2 Complete Cox · eTable3 Stratified E-values · eTable4 Model Specification · eTable5 Leave-one-out · eTable6 Geographic Variation · eTable7 Schoenfeld + Landmark · eTable8 Hyperparameters · eTable9 TARGET Protocol Specification · eTable10 Decision Matrix · eTable11 Included vs Excluded

---

**Data Sharing Statement**: Data used in this study are publicly available from SEER (seer.cancer.gov, Research Plus Data, November 2024 submission), TCGA-LIHC (portal.gdc.cancer.gov), and ICGC-LIRI-JP (dcc.icgc.org). Analysis code is available at github.com/woodhaha/hepatobiliary-cancer-seer-analysis.

**Author Contributions**: Drs Xu and Bai contributed equally as co-first authors. Dr Zhou had full access to all data and takes responsibility for data integrity and analysis. Study concept and design: Xu, Zhou. Data acquisition: Xu, Bai, Cai. Statistical analysis: Zhou. Machine learning: Zhou, Bai. Manuscript drafting: Zhou. Critical revision: all authors. Study supervision: Zhou.

**Conflict of Interest Disclosures**: None reported.

**Funding/Support**: This research did not receive any specific grant from funding agencies in the public, commercial, or not-for-profit sectors.

**Role of the Funder/Sponsor**: Not applicable (no funding).

**Additional Contributions**: The authors acknowledge the National Cancer Institute for the SEER Program, the TCGA Research Network, and the ICGC consortium for making data publicly available.

**Meeting Presentation**: Not presented.
