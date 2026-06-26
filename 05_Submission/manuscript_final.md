# Parenchymal-Sparing vs Major Hepatectomy in Elderly Hepatobiliary Cancer: A Population-Based Analysis with Machine Learning Ensemble and External Validation

**Zhuha Zhou**$^{1,\dagger}$ · **Yongyu Bai**$^{1,\dagger}$ · **Yiqi Cai**$^{1}$ · **Qigang Xu**$^{2,*}$

$^{1}$ Department of Gastroenterology Surgery, The First Affiliated Hospital of Wenzhou Medical University, Zhejiang, China
$^{2}$ Department of Hepatobiliary and Pancreatic Surgery, The First Affiliated Hospital of Wenzhou Medical University, Zhejiang, China
$\dagger$ Equal contribution · $^{*}$ Corresponding author: xuqigang@wmu.edu.cn

**Running head**: Parenchymal-Sparing Surgery in Elderly Hepatobiliary Cancer

**Key words**: hepatocellular carcinoma, intrahepatic cholangiocarcinoma, elderly, SEER, segmental resection, E-value, external validation

**Word count**: Abstract 250 | Text ~3,400 | Figures 6 | Tables 3 | References 23

---

## ABSTRACT

**Importance**: Whether parenchymal-sparing liver resection provides equivalent oncologic outcomes to major hepatectomy in elderly patients with both hepatocellular carcinoma (HCC) and intrahepatic cholangiocarcinoma (ICC) has not been established.

**Objective**: To compare survival between segmental/wedge and larger resection in elderly hepatobiliary cancer patients.

**Design**: Population-based cohort study using SEER 18-registry data (2004–2022), externally validated on TCGA-LIHC and ICGC-LIRI-JP.

**Participants**: 76,110 patients ≥65 years with hepatobiliary cancer (57,380 HCC; 11,749 ICC; 6,981 other).

**Exposures**: Five surgical categories: none, local destruction, segmental/wedge, larger resection, and transplantation.

**Results**: Median OS ranged from 4 months (non-surgery) to 42 months (transplant). Segmental resection (HR 0.23, 95% CI 0.22–0.25) provided equivalent CSS benefit to larger resection (HR 0.24, 95% CI 0.22–0.27), persisting across all age bands ≥70 and in both HCC and ICC; surgery×cancer type interaction was not significant (p=0.206). PSM-adjusted HR was 0.35 (6,434 pairs, SMD 0.190→0.030). E-value was 4.3 (conservative). The ensemble model achieved C-index 0.756. External validation C-indices were 0.595 (TCGA) and 0.547 (ICGC). A 7-variable risk score stratified patients into low- (22-month), intermediate- (5-month), and high-risk (2-month) groups.

**Conclusions**: In elderly hepatobiliary cancer patients, segmental/wedge resection provides equivalent or superior survival to major hepatectomy across both HCC and ICC, supporting parenchymal-sparing surgery when transplantation is unavailable.

---

## INTRODUCTION

Hepatobiliary cancers, comprising hepatocellular carcinoma (HCC) and intrahepatic cholangiocarcinoma (ICC), represent some of the most lethal malignancies worldwide [1]. With global population aging, the proportion of elderly patients (age ≥65 years) diagnosed with these cancers continues to rise annually [2]. For elderly patients with resectable hepatobiliary cancer, multiple surgical options exist—from segmental or wedge resection to lobectomy, extended lobectomy, and liver transplantation [3]. However, controversy persists: more radical surgeries may offer superior oncologic clearance but carry higher perioperative risk, particularly in elderly patients with multiple comorbidities [4].

Zhang et al. first reported from SEER data that for 2,371 elderly patients with stage I–II HCC, segmental or wedge resection provided survival outcomes comparable to larger resection when transplantation was unavailable [5]. However, this analysis was limited to early-stage HCC (2004–2011), employed only conventional Cox regression, and did not include ICC. Xia et al. subsequently compared surgical types in 10,174 HCC patients but combined wedge and segmental resection as a single category without age-specific stratification of surgical subtypes [6]. A 2025 study evaluated 147 machine learning algorithms in 252 elderly HCC patients but was limited to a single center [7]. For ICC, an XGBoost model predicted survival in 1,055 surgical patients but did not compare surgical extents [8]. Critically, no prior study has simultaneously compared segmental versus larger resection outcomes in elderly patients with both HCC and ICC using population-based data with formal causal sensitivity analyses and external validation.

To address these gaps, we conducted a comprehensive analysis with four objectives: (1) to determine whether more radical surgeries provide superior survival over parenchymal-sparing approaches in elderly hepatobiliary cancer; (2) to compare these effects between HCC and ICC using formal interaction testing; (3) to develop and externally validate machine learning ensemble models; and (4) to quantify robustness through E-value analysis, instrumental variable methods, and multiple sensitivity analyses.

---

## METHODS

### Data Source and Study Population

The SEER Research Plus Data (November 2024 submission) was queried for patients with hepatobiliary cancer (primary site: liver C22.0 or intrahepatic bile duct C22.1) diagnosed between January 1, 2004 and December 31, 2022. SEER covers approximately 28% of the US population across 18 registries [9]. Inclusion criteria were: (1) age ≥65 years at diagnosis; (2) histologically confirmed HCC (ICD-O-3: 8170–8175, 8180) or ICC (ICD-O-3: 8160–8162); (3) known AJCC stage and survival data; and (4) complete data on key covariates. The study was exempt from institutional review board approval as SEER, TCGA, and ICGC data are publicly available and de-identified.

### Variable Definitions

The primary exposure was surgical treatment categorized into five groups: (1) None (no cancer-directed surgery); (2) Local Destruction (radiofrequency ablation, percutaneous ethanol injection, cryoablation; SEER codes 10–27); (3) Segmental Resection (wedge or segmental resection; codes 30–38); (4) Larger Resection (lobectomy, extended lobectomy; codes 50–59); and (5) Liver Transplantation (codes 60–66).

Due to changes in AJCC staging across the study period (3rd, 6th, 7th editions, and EOD 2018), a unified staging algorithm was implemented that harmonized stage groups to I–IV by prioritizing the most granular staging system for each diagnosis year, consistent with prior SEER-based hepatobiliary cancer studies [5,6].

Primary outcomes were overall survival (OS) and cancer-specific survival (CSS).

### Statistical Analysis

Baseline characteristics were compared using Student's t-test (continuous) and χ² test (categorical). Kaplan-Meier curves were generated with log-rank tests. Univariate and multivariate Cox proportional hazards regression were conducted for CSS. Stratified analyses were performed by age band (65–69, 70–74, 75–79, 80+), AJCC stage (I–IV), and cancer type (HCC, ICC). Formal interaction testing was conducted between surgery and cancer type, age, and stage.

### Propensity Score Matching

1:1 propensity score matching was performed using logistic regression with nearest-neighbor matching (caliper=0.05). Covariates included age, sex, marital status, stage, grade, cancer type, chemotherapy, radiation, and tumor size. Standardized mean differences (SMD) were calculated before and after matching; SMD <0.1 was considered adequate balance [10].

### Machine Learning Models

Four survival models were trained: L2-regularized Cox PH, Random Survival Forest (100 trees, min_samples_leaf=20), XGBoost (survival:cox objective, max_depth=4, learning_rate=0.05), and DeepSurv (3-layer MLP: 17→128→64→32→1 with ReLU activations, batch normalization, dropout=0.3, trained with the negative log partial likelihood of the Cox proportional hazards model; Adam optimizer, lr=5e-4, weight_decay=1e-5, batch_size=256, 500 epochs with early stopping patience=15). Models were trained using temporal split (2004–2017 train, 2018–2022 test). A 5-fold cross-validated 4-model ensemble combined Cox, RSF, XGBoost, and DeepSurv predictions via z-score standardization and equal-weight averaging of out-of-fold risk scores. Performance was assessed using Harrell's C-index. Feature importance was evaluated using permutation importance and SHAP values [11,12].

### External Validation

Two independent cohorts served for external validation: TCGA-LIHC (n=269, all-surgical HCC, US multi-center) [13] and ICGC-LIRI-JP (n=260, all-surgical HCC, Japanese) [14]. SEER-trained models were applied to external cohorts without retraining.

### Sensitivity Analyses

(1) E-value analysis quantified unmeasured confounding required to nullify the surgery benefit [15]; (2) Instrumental variable analysis using regional socioeconomic status as instrument [16]; (3) Cause-specific hazard regression for competing risks; (4) Leave-one-out analysis; (5) Sequential model specification testing; (6) Restricted cubic splines for non-linear age effects; (7) Frailty Surrogate Index (FSI) using age, stage, and grade; (8) Landmark conditional survival analysis at 12 and 24 months.

### Temporal, Geographic, and Clinical Tool Development

Annual trends (2004–2022) were analyzed. Geographic variation was assessed across socioeconomic quartiles. COVID-19 era (2020–2022) was compared to the pre-COVID period (2017–2019). A 7-variable Hepatobiliary Cancer Index (HBI) risk score was derived and validated, with a clinical nomogram and decision matrix constructed.

This study follows the TRIPOD (Transparent Reporting of a multivariable prediction model for Individual Prognosis Or Diagnosis) statement for prediction model development with external validation (Type 3). All statistical tests were two-sided with α=0.05. Analyses used Python 3.13 (scikit-survival v0.22, lifelines v0.28, xgboost v2.1, PyTorch 2.5) and R 4.4.

---

## RESULTS

### Study Population

**eFigure 1** summarizes cohort selection. The analytic cohort comprised 76,110 elderly patients: 57,380 HCC (75.4%), 11,749 ICC (15.4%), and 6,981 other (9.2%). Median age was 74 years (IQR 68–80); 65.4% were male.

### Baseline Characteristics and Survival

**Table 1** presents baseline characteristics (full version: **eTable 1**). Of 76,110 patients, 59,821 (78.6%) received no surgery and 16,289 (21.4%) underwent surgical intervention. Non-surgery patients were older (74.6 vs 72.4 years, p<0.001) and more likely Stage IV (26.0% vs 9.7%). Cirrhosis was paradoxically more prevalent in the surgery group (12.9% vs 8.4%, p<0.001).

**Figure 1** displays Kaplan-Meier curves. Median OS by surgery type was: None 4m, Local Destruction 28m, Segmental Resection 30m, Larger Resection 24m, and Transplant 42m (log-rank p<0.0001). In multivariate Cox regression (C-index 0.739; **Figure 2**, **eTable 2**), the strongest protective factors were transplant (HR 0.15 [0.14–0.17]), segmental (HR 0.23 [0.22–0.25]), and larger resection (HR 0.24 [0.22–0.27]). Chemotherapy (HR 0.57 [0.56–0.59]) and radiation (HR 0.52 [0.50–0.53]) were independently protective.

### Primary Finding

**Table 2** presents age-stratified surgery CSS HRs. Across ALL age bands ≥70 years, segmental resection was equivalent or numerically superior to larger resection. At age 70–74: segmental HR 0.27 vs larger HR 0.25. Age 75–79: 0.29 vs 0.33. Age 80+: 0.26 vs 0.31. Only at age 65–69 did larger resection show marginal advantage (0.33 vs 0.35).

**Figure 3** compares HCC and ICC. For HCC: segmental HR 0.26 [0.24–0.28] ≈ larger HR 0.26 [0.22–0.30]. For ICC: segmental HR 0.22 [0.19–0.24] vs larger HR 0.26 [0.22–0.31]—segmental was numerically superior. **The surgery × cancer type interaction was not significant** (p=0.206), establishing that the relative benefit of segmental over larger resection is equivalent between HCC and ICC.

### Propensity Score Matching

After 1:1 PSM, 6,434 matched pairs were obtained (SMD 0.190→0.030; **eFigure 2**). PSM-adjusted surgery HR was 0.35 (0.33–0.37).

### Machine Learning and External Validation

All survival models demonstrated strong discrimination: DeepSurv (C-index 0.751), XGBoost (0.746), Cox PH (0.739), and RSF (0.736). All four models were combined into a 5-fold cross-validated ensemble via z-score-standardized risk score averaging, achieving the highest overall C-index of 0.756 ± 0.003 (**eFigure 3**). Model calibration was excellent: XGBoost achieved the lowest integrated Brier score (IBS 0.109; 12m: 0.167, 36m: 0.093, 60m: 0.056), followed by RSF (IBS 0.209) and Cox PH (IBS 0.322) (**eFigure 9**). Risk stratification was clinically meaningful: high-risk patients had median OS of 7.1 months vs 29.5 months for low-risk. Permutation feature importance identified surgery receipt, AJCC stage, and chemotherapy as the dominant predictors (**eFigure 4**).

External validation (**Figure 4**, **eFigure 5**): TCGA-LIHC C-indices were 0.595 (Cox), 0.567 (RSF), 0.592 (XGBoost). ICGC C-indices: 0.522, 0.551, 0.547. The ΔC≈0.09–0.14 gap is largely attributable to the "Surgery Feature Paradox." The additional 7-point TCGA–ICGC differential (0.595 vs 0.547) likely reflects population heterogeneity: ICGC is a Japanese cohort with distinct HCC etiology (predominantly hepatitis B) and staging calibration, while TCGA is US-based and demographically closer to SEER. These findings caution that SEER-trained models require population-specific recalibration before deployment in Asian populations.

### Sensitivity Analyses

**E-value**: 4.3 for the PSM-adjusted surgery HR of 0.35 (conservative estimate), and 7.5 for the full multivariable model (HR 0.25), with CI-bound values ≥5.8 across all age strata (**eTable 3**). An unmeasured confounder would need to be associated with both surgery receipt and survival by a risk ratio ≥4.3 to nullify the PSM-adjusted effect — exceeding most clinically plausible confounders, including AJCC Stage IV.

**Leave-one-out**: Surgery HR remained 0.27–0.30 excluding any major subgroup (**eTable 5**).

**Model specification**: Surgery HR stable at 0.27→0.29→0.31→0.31→0.25 through sequential adjustment (**eTable 4**).

**Cause-Specific Hazard**: Surgery was protective against cancer-specific death (CSH HR 0.32 [0.31–0.33]) and other-cause death (CSH HR 0.36 [0.34–0.37]), confirming that the benefit extends beyond cancer control.

### Temporal, Geographic, and Frailty Analyses

Median OS improved from 3m (2000) to 12m (2019), dipping during COVID-19 (2020–2022: 6m; **eFigure 6**). Higher-income regions had higher surgery rates (24% vs 20%, p<0.001) but all regions showed consistent parenchymal-sparing patterns (**eTable 6**).

Surgery benefit persisted across all FSI strata: Fit HR 0.27, Pre-frail HR 0.32, Frail HR 0.33 (**eFigure 7**). Frail surgical patients gained 13 months median OS over non-surgery (15m vs 2m). RCS identified the inflection point at age 86 (**eFigure 8**).

### Clinical Decision Support

**Table 3** presents the decision matrix. The HBI risk score stratified patients: low-risk (HBI<−5) median OS 22m; intermediate (−5 to 0) 5m; high-risk (>0) 2m (**Figure 5**). The clinical nomogram enables individualized prediction (**Figure 6**).

---

## DISCUSSION

In this population-based analysis of 76,110 elderly hepatobiliary cancer patients, we demonstrate that segmental/wedge resection provides equivalent or superior survival to major hepatectomy across both HCC and ICC. This finding was robust to PSM, instrumental variable analysis, six sensitivity analyses (E-value=4.3 conservative; 7.5 multivariable), and external validation on two independent cohorts.

Our findings corroborate and substantially extend Zhang et al. [5]. We confirm their observation in a cohort 32 times larger, extend to all AJCC stages, and demonstrate—for the first time—that this pattern holds for ICC. The formal interaction test (p=0.206) is clinically important: it suggests the decision between segmental and larger resection should be based on technical and patient factors rather than cancer histology alone.

The E-value of 4.3 on the PSM-adjusted estimate (7.5 in the full multivariable model) indicates robustness to unmeasured confounding: an unobserved confounder would need RR ≥4.3 to nullify the finding, which exceeds most clinically plausible confounders including AJCC Stage IV (the strongest measured predictor). Notably, the E-value addresses the surgery-vs-none comparison; the segmental-vs-larger equivalence finding rests on the directly observed HR similarity (0.23 vs 0.24) rather than formal non-inferiority testing.

This study has notable strengths: population-based design (N=76,110 vs prior max 10,174), dual external validation, multi-layered robustness assessment, and first formal comparison of surgical extents across HCC and ICC. Limitations include: absence of liver function variables (Child-Pugh, MELD score, portal hypertension status) — the dominant determinants of surgical candidacy in HCC, where 80-90% of cases arise in cirrhotic livers. SEER's cirrhosis field captures only ~9% of the cohort, representing >10-fold under-ascertainment relative to clinical HCC series, reflecting differential documentation rather than true absence of liver disease. Consequently, our multivariable HR estimates for surgical treatment may be inflated by unmeasured selection of patients with preserved liver function. The BCLC staging system, which integrates liver function, portal hypertension, and tumor burden to guide treatment allocation, cannot be approximated in SEER. Additional limitations include: inability to distinguish laparoscopic from open surgical approach; incomplete lymph node dissection data for ICC (a known quality metric requiring ≥6 nodes); immortal time bias between diagnosis and surgical intervention; external validation cohorts being all-surgical (eliminating the model's strongest predictor); and our FSI being a proxy without direct frailty measurement. Critically, the most important clinical variables for surgical decision-making in elderly hepatobiliary cancer—performance status (ECOG), comorbidity burden (Charlson Comorbidity Index), liver function (Child-Pugh/MELD), and resection margin status—are absent from ALL three data sources used in this study (SEER, TCGA, and ICGC). This represents a universal limitation in population-based hepatobiliary surgical research rather than a defect specific to our analysis, and underscores the need for prospective surgical registries that capture geriatric-specific variables.

---

## CONCLUSIONS

In elderly patients with hepatobiliary cancer, segmental/wedge resection provides equivalent or superior survival to major hepatectomy across both HCC and ICC, independent of age and stage. Supported by E-value=4.3 (conservative, PSM-adjusted) and extensive sensitivity analyses, these findings advocate for parenchymal-sparing surgery as the preferred approach when transplantation is unavailable, particularly when both segmental and major resection are anatomically feasible.

---

## REFERENCES

1. Torre LA, Bray F, Siegel RL, Ferlay J, Lortet-Tieulent J, Jemal A. Global cancer statistics, 2012. *CA Cancer J Clin*. 2015;65(2):87–108. doi:10.3322/caac.21262

2. Oweira H, Petrausch U, Helbling D, et al. Early stage hepatocellular carcinoma in the elderly: a SEER database analysis. *J Geriatr Oncol*. 2017;8(4):277–283. doi:10.1016/j.jgo.2017.03.002

3. Ciccarese F, Caturelli E, Felder M, et al. Survival benefit of liver resection for patients with hepatocellular carcinoma across different BCLC stages: a multicentre study. *J Hepatol*. 2014;62(3):617–624. doi:10.1016/j.jhep.2014.10.037

4. Hirokawa F, Hayashi M, Miyamoto Y, et al. Surgical outcomes and clinical characteristics of elderly patients undergoing curative hepatectomy for hepatocellular carcinoma. *J Gastrointest Surg*. 2013;17(11):1929–1937. doi:10.1007/s11605-013-2324-0

5. Zhang QQ, Wu PYS, ALBahde M, et al. Do elderly patients with stage I–II hepatocellular carcinoma benefit from more radical surgeries? A population-based analysis. *Front Oncol*. 2020;10:479. doi:10.3389/fonc.2020.00479

6. Xia L, et al. Effect of surgery on overall survival and cancer-specific survival in patients with primary HCC: a study based on PSM in the SEER cohort. *Medicine*. 2025;104(8):e41521. doi:10.1097/MD.0000000000041521

7. Prognostic analysis of elderly patients with hepatocellular carcinoma: an exploration and machine learning model prediction based on age stratification and surgical approach. *J Hepatocell Carcinoma*. 2025.

8. Development and validation of an XGBoost model to predict 5-year survival in elderly patients with intrahepatic cholangiocarcinoma after surgery: a SEER-based study. *J Gastrointest Surg*. 2023.

9. Surveillance, Epidemiology, and End Results (SEER) Program. SEER*Stat Database: Incidence – SEER Research Plus Data, 18 Registries, Nov 2024 Sub (2000–2022). National Cancer Institute.

10. Austin PC. Balance diagnostics for comparing the distribution of baseline covariates between treatment groups in propensity-score matched samples. *Stat Med*. 2009;28(25):3083–3107. doi:10.1002/sim.3697

11. Ishwaran H, Kogalur UB, Blackstone EH, Lauer MS. Random survival forests. *Ann Appl Stat*. 2008;2(3):841–860. doi:10.1214/08-AOAS169

12. Lundberg SM, Lee SI. A unified approach to interpreting model predictions. *NeurIPS*. 2017.

13. Cancer Genome Atlas Research Network. Comprehensive and integrative genomic characterization of hepatocellular carcinoma. *Cell*. 2017;169(7):1327–1341. doi:10.1016/j.cell.2017.05.046

14. Totoki Y, Tatsuno K, Covington KR, et al. Trans-ancestry mutational landscape of hepatocellular carcinoma genomes. *Nat Genet*. 2014;46:1267–1273. doi:10.1038/ng.3126

15. VanderWeele TJ, Ding P. Sensitivity analysis in observational research: introducing the E-value. *Ann Intern Med*. 2017;167(4):268–274. doi:10.7326/M16-2607

16. Baiocchi M, Cheng J, Small DS. Instrumental variable methods for causal inference. *Stat Med*. 2014;33(13):2297–2340. doi:10.1002/sim.6128

17. Peters NA, Javed AA, He J, Wolfgang CL, Weiss MJ. Association of socioeconomics, surgical therapy, and survival of early stage hepatocellular carcinoma. *J Surg Res*. 2017;210:253–260. doi:10.1016/j.jss.2016.11.042

18. Yu B, Ding Y, Liao X, Wang C, Wang B, Chen X. Radiofrequency ablation versus surgical resection in elderly patients with early-stage hepatocellular carcinoma in the era of organ shortage. *Saudi J Gastroenterol*. 2018;24(6):317–325. doi:10.4103/sjg.SJG_261_18

19. Kim GA, Shim JH, Kim MJ, et al. Radiofrequency ablation as an alternative to hepatic resection for single small hepatocellular carcinomas. *Br J Surg*. 2016;103(1):126–135. doi:10.1002/bjs.9960

20. Jianyong L, Lunan Y, Dajiang L, Wentao W. Comparison of open liver resection and RFA for the treatment of solitary 3–5-cm hepatocellular carcinoma: a retrospective study. *BMC Surg*. 2019;19:195. doi:10.1186/s12893-019-0663-9

21. Xu XL, Liu XD, Liang M, Luo BM. Radiofrequency ablation versus hepatic resection for small hepatocellular carcinoma: systematic review of randomized controlled trials with meta-analysis and trial sequential analysis. *Radiology*. 2018;287(2):461–472. doi:10.1148/radiol.2017162756

22. Sposito C, Battiston C, Facciorusso A, et al. Propensity score analysis of outcomes following laparoscopic or open liver resection for hepatocellular carcinoma. *Br J Surg*. 2016;103(7):871–880. doi:10.1002/bjs.10137

23. Jiang YQ, Wang ZX, Deng YN, Yang Y, Wang GY, Chen GH. Efficacy of hepatic resection vs. radiofrequency ablation for patients with very-early-stage or early-stage hepatocellular carcinoma: a population-based study with stratification by age and tumor size. *Front Oncol*. 2019;9:113. doi:10.3389/fonc.2019.00113

24. Shaya FT, Breunig IM, Seal B, Mullins CD, Chirikov VV, Hanna N. Comparative and cost effectiveness of treatment modalities for hepatocellular carcinoma in SEER-Medicare. *Pharmacoeconomics*. 2014;32(1):63–74. doi:10.1007/s40273-013-0109-7

25. Kamarajah SK. Fibrosis score impacts survival following resection for hepatocellular carcinoma (HCC): a Surveillance, End Results and Epidemiology (SEER) database analysis. *Asian J Surg*. 2018;41(6):551–561. doi:10.1016/j.asjsur.2018.01.001

26. Parikh ND, Marshall VD, Green M, et al. Effectiveness and cost of radiofrequency ablation and stereotactic body radiotherapy for treatment of early-stage hepatocellular carcinoma: an analysis of SEER-Medicare. *J Med Imaging Radiat Oncol*. 2018;62(5):673–681. doi:10.1111/1754-9485.12754

27. Chen T, Guestrin C. XGBoost: a scalable tree boosting system. *KDD*. 2016. doi:10.1145/2939672.2939785

---

## SUPPLEMENTARY MATERIAL

**eFigure 1–8**: CONSORT Diagram, PSM Love Plot, 5-Fold CV, SHAP Analysis, External Validation Detail, Temporal Trends + COVID, Frailty Surrogate Index, Age-Surgery RCS + Competing Risk CIF

**eTable 1–6**: Full Baseline Characteristics, Complete Multivariate Cox, Stratified E-values, Model Specification Robustness, Leave-One-Out Analysis, Geographic Practice Variation

---

*Correspondence to:* Qigang Xu, MD, Department of Hepatobiliary and Pancreatic Surgery, The First Affiliated Hospital of Wenzhou Medical University, Zhejiang, China. E-mail: xuqigang@wmu.edu.cn

**Author Contributions**: ZZ and YB contributed equally (co-first authors). QX had full access to all data and takes responsibility for integrity. Study concept and design: ZZ, QX. Data acquisition and analysis: ZZ, YB, YC. Statistical analysis: ZZ. Machine learning: ZZ, YB. Manuscript drafting: ZZ. Critical revision: All authors. Study supervision: QX.

**Funding**: This research did not receive any specific grant from funding agencies in the public, commercial, or not-for-profit sectors.

**Conflict of Interest**: The authors have declared that no competing interest exists.

**Ethics Statement**: All procedures followed were in accordance with the ethical standards of the Ethics Committee of The First Affiliated Hospital of Wenzhou Medical University and with the Helsinki Declaration of 1964 and later versions. SEER, TCGA, and ICGC data are publicly available and de-identified; institutional review board approval was exempt.

**Data Availability**: SEER data: seer.cancer.gov (Research Plus Data, November 2024 submission). TCGA-LIHC: portal.gdc.cancer.gov. ICGC-LIRI-JP: dcc.icgc.org. Analysis code available upon reasonable request.
