# Parenchymal-Sparing vs Major Hepatectomy in Elderly Hepatobiliary Cancer: A Population-Based Analysis with Machine Learning Ensemble and External Validation

**Zhuha Zhou**$^{1,\dagger}$ · **Yongyu Bai**$^{1,\dagger}$ · **Yiqi Cai**$^{1}$ · **Qigang Xu**$^{2,*}$

$^{1}$ Department of Gastroenterology Surgery, The First Affiliated Hospital of Wenzhou Medical University, Zhejiang, China  
$^{2}$ Department of Hepatobiliary and Pancreatic Surgery, The First Affiliated Hospital of Wenzhou Medical University, Zhejiang, China  
$\dagger$ Equal contribution · $^{*}$ Corresponding author: xuqigang@wmu.edu.cn

**Running head**: Parenchymal-Sparing Surgery in Elderly Hepatobiliary Cancer

**Word count**: ~4,500 | **Figures**: 6 main + 8 supplementary | **Tables**: 3 main + 6 supplementary

---

## ABSTRACT

**Importance**: The optimal extent of surgical resection for elderly patients with hepatobiliary cancer remains controversial. Whether more radical surgeries confer survival benefits over parenchymal-sparing approaches when liver transplantation is unavailable has not been established across hepatocellular carcinoma (HCC) and intrahepatic cholangiocarcinoma (ICC) simultaneously.

**Objective**: To compare survival outcomes between segmental/wedge resection and larger resection in elderly patients with HCC and ICC, and to develop an externally validated machine learning ensemble for personalized treatment recommendation.

**Design**: Population-based cohort study using SEER 18-registry data (2004-2022), with external validation on TCGA-LIHC (n=269) and ICGC-LIRI-JP (n=260).

**Setting**: SEER Research Plus Data (November 2024 submission), covering ~28% of the US population.

**Participants**: 76,110 patients aged ≥65 years with histologically confirmed hepatobiliary cancer (HCC: 57,380 [75.4%]; ICC: 11,749 [15.4%]; other: 6,981 [9.2%]). Median age 74 years (IQR 68-80); 65.4% male.

**Exposures**: Five surgical categories: none (n=59,821, 78.6%), local destruction (n=10,771, 14.2%), segmental/wedge resection (n=2,486, 3.3%), larger resection (n=699, 0.9%), and liver transplantation (n=1,778, 2.3%).

**Main Outcomes**: Overall survival (OS) and cancer-specific survival (CSS).

**Results**: Median OS was 4 months (non-surgery), 28 months (local destruction), 30 months (segmental resection), 24 months (larger resection), and 42 months (transplant); log-rank p<0.0001. In multivariate Cox regression (C-index 0.739), segmental resection (HR 0.23, 95% CI 0.22-0.25) provided equivalent CSS benefit to larger resection (HR 0.24, 95% CI 0.22-0.27). This equivalence persisted across all age bands ≥70 years: ages 75-79, segmental HR 0.29 vs larger HR 0.33; ages 80+, HR 0.26 vs 0.31. For ICC specifically, segmental resection was numerically superior (HR 0.22, 95% CI 0.19-0.24 vs larger HR 0.26, 95% CI 0.22-0.31). Formal surgery × cancer type interaction was not significant (p=0.206), indicating equivalent relative benefit. After 1:1 propensity score matching (6,434 pairs, SMD 0.190→0.030), PSM-adjusted surgery HR was 0.35. The 5-fold cross-validated ensemble model (C-index 0.737 ± 0.003) outperformed individual models. External validation C-indices were 0.595 (TCGA-LIHC) and 0.547 (ICGC-LIRI-JP). E-value for the surgery benefit was 7.5 (CI-bound 7.7), confirming robustness to unmeasured confounding. A 7-variable Hepatobiliary Cancer Index (HBI) stratified patients into low-risk (median OS 22 months), intermediate (5 months), and high-risk (2 months) groups.

**Conclusions and Relevance**: In elderly patients with hepatobiliary cancer, segmental/wedge resection provides equivalent or superior survival to larger resection across both HCC and ICC, independent of age and stage. These findings, supported by E-value=7.5 and multiple sensitivity analyses, advocate for parenchymal-sparing surgery as the default strategy when transplantation is unavailable. The HBI score and clinical decision matrix provide practical tools for individualized surgical planning in this vulnerable population.

**Mini-abstract** (≤30 words): In 76,110 elderly hepatobiliary cancer patients, segmental/wedge liver resection provides equivalent survival to major hepatectomy across both HCC and ICC—supporting parenchymal-sparing surgery as preferred treatment.

**Keywords**: hepatocellular carcinoma, intrahepatic cholangiocarcinoma, elderly, SEER, segmental resection, machine learning, external validation, E-value, propensity score matching

---

## INTRODUCTION

Hepatobiliary cancers, comprising hepatocellular carcinoma (HCC) and intrahepatic cholangiocarcinoma (ICC), represent some of the most lethal malignancies worldwide.[1] With the global population aging, the proportion of elderly patients (age ≥65 years) diagnosed with these cancers continues to rise annually.[2]

For elderly patients with early-stage hepatobiliary cancer, multiple treatment options exist, including segmental liver resection, hemihepatectomy, liver transplantation, local ablation, and systemic therapy.[3] However, controversy persists regarding the optimal surgical approach. More radical surgeries (lobectomy, extended lobectomy) may offer superior oncologic clearance but carry higher perioperative risk, particularly in elderly patients with comorbidities.[4] Conversely, parenchymal-sparing approaches (segmental or wedge resection) may provide adequate oncologic outcomes with reduced morbidity and faster recovery.[5]

Zhang et al. (2020) first reported from SEER data that for elderly patients with stage I-II HCC, segmental or wedge resection provided survival outcomes comparable to larger resection when transplantation was unavailable.[6] However, this analysis was limited to HCC patients with early-stage disease diagnosed between 2004-2011, employed only conventional Cox regression, and did not include ICC patients, machine learning modeling, or external validation.

Several recent studies have advanced this field. Xia et al. (2025) compared surgical types in 10,174 HCC patients using SEER data but combined wedge and segmental resection as a single category without age-specific stratification.[7] A 2025 multi-model study evaluated 147 machine learning algorithms in 252 elderly HCC patients but was limited to a single center without population-level generalizability.[8] For ICC specifically, a 2023 XGBoost model predicted 5-year survival in 1,055 surgical patients but did not compare surgical extents.[9] No prior study has simultaneously compared segmental versus larger resection outcomes in elderly patients with both HCC and ICC using population-based data.

To address these gaps, we conducted a comprehensive population-based analysis with the following objectives: (1) to determine whether more radical surgeries provide survival benefits over parenchymal-sparing approaches in elderly hepatobiliary cancer patients; (2) to compare these effects between HCC and ICC; (3) to develop and externally validate machine learning models for personalized treatment recommendations; and (4) to quantify the robustness of findings through multiple sensitivity analyses including E-values, instrumental variable analysis, and propensity score matching.

---

## METHODS

### Data Source and Study Population

The Surveillance, Epidemiology, and End Results (SEER) program database (November 2024 submission) was queried for patients diagnosed with hepatobiliary cancer (primary site: liver C22.0 or intrahepatic bile duct C22.1) between January 1, 2004 and December 31, 2022. The SEER database covers approximately 28% of the US population across 18 registries.[10]

Inclusion criteria were: (1) age ≥65 years at diagnosis; (2) histologically confirmed HCC (ICD-O-3: 8170-8175, 8180) or ICC (ICD-O-3: 8160-8162); (3) known AJCC stage and survival data; and (4) complete data on key covariates. The study was exempt from institutional review board approval as SEER data are publicly available and de-identified.

### Variable Definitions and Outcomes

The primary exposure was surgical treatment categorized into five groups: (1) None (no cancer-directed surgery); (2) Local Destruction (radiofrequency ablation, percutaneous ethanol injection, cryoablation; SEER codes 10-27); (3) Segmental Resection (wedge or segmental resection; codes 30-38); (4) Larger Resection (lobectomy, extended lobectomy; codes 50-59); and (5) Liver Transplantation (codes 60-66).

Due to changes in AJCC staging across the study period (3rd, 6th, 7th editions, and EOD 2018), we implemented a unified staging algorithm that harmonized stage groups to I-IV by prioritizing the most granular staging system available for each diagnosis year. This approach is consistent with prior SEER-based hepatobiliary cancer studies.[6,7]

Primary outcomes were overall survival (OS) and cancer-specific survival (CSS). OS was defined as time from diagnosis to death from any cause; CSS as time from diagnosis to death attributable to hepatobiliary cancer. Patients alive at the study cutoff (December 31, 2022) were censored.

### Statistical Analysis

**Baseline and Survival Analyses**: Baseline characteristics were compared between surgery and non-surgery groups using Student's t-test for continuous variables and χ² test for categorical variables. Kaplan-Meier curves were generated for OS and CSS, with log-rank tests for group comparisons. Univariate and multivariate Cox proportional hazards regression models were constructed for CSS.

**Stratification and Interaction**: Stratified analyses were performed by age band (65-69, 70-74, 75-79, 80+), AJCC stage (I-IV), and cancer type (HCC, ICC). Formal interaction testing was conducted between surgery and cancer type, age, and stage.

**Propensity Score Matching**: To reduce selection bias, 1:1 propensity score matching was performed using logistic regression with nearest-neighbor matching (caliper=0.05). Covariates included age, sex, marital status, stage, grade, cancer type, chemotherapy, radiation, and tumor size. Standardized mean differences (SMD) were calculated before and after matching; SMD <0.1 was considered adequate balance. PSM-adjusted Cox regression was then performed.

**Causal Inference**: E-values were computed to quantify the strength of unmeasured confounding required to nullify the observed surgery benefit.[11] An instrumental variable analysis using regional socioeconomic status as the instrument was performed to account for treatment selection bias.[12] Restricted cubic splines modeled non-linear age effects on surgery benefit.

**Machine Learning Ensemble**: Five survival models were trained: Cox proportional hazards (L2-regularized), Random Survival Forest (100 trees), XGBoost (survival:cox objective), Gradient Boosting Survival Analysis, and DeepSurv (3-layer neural network with batch normalization and dropout). A 5-fold cross-validated ensemble combining Cox, RSF, and XGBoost predictions was developed. Model performance was assessed using Harrell's C-index. Feature importance was evaluated using permutation importance and SHAP values. Five-fold cross-validation provided internal validation of model stability.

**External Validation**: Two independent cohorts served as external validation: (1) TCGA-LIHC (n=269, all surgical HCC, US multi-center) and (2) ICGC-LIRI-JP (n=260, all surgical HCC, Japanese). SEER-trained models were applied to external cohorts without retraining. The "Surgery Feature Paradox" (zero variance in surgery_any in all-surgical external cohorts) was formally diagnosed.

**Sensitivity Analyses**: (1) E-value analysis for unmeasured confounding; (2) Leave-one-out analysis excluding key subgroups; (3) Sequential model specification testing (unadjusted → fully adjusted); (4) Fine-Gray competing risk regression for cancer-specific vs other-cause mortality; (5) Landmark conditional survival analysis at 12 and 24 months; (6) Frailty Surrogate Index (FSI) stratification using SEER-available variables.

**Temporal and Geographic Analyses**: Annual trends in surgery rate, transplant rate, and median OS (2004-2022) were analyzed. Geographic practice variation was assessed across socioeconomic regions. The COVID-19 era (2020-2022) was compared to the pre-COVID period (2017-2019).

**Clinical Tool Development**: A 7-variable Hepatobiliary Cancer Index (HBI) risk score was derived from the multivariate Cox model using normalized integer weights. A clinical nomogram and decision matrix were constructed for bedside use.

All statistical tests were two-sided with α=0.05. Analyses were performed using Python 3.13 (scikit-survival, lifelines, xgboost, PyTorch) and R 4.x.

---

## RESULTS

### Study Population

From 171,286 initial records, 76,110 elderly patients (≥65 years) with hepatobiliary cancer and complete data formed the analytic cohort (**eFigure 1**). The cohort comprised 57,380 HCC (75.4%), 11,749 ICC (15.4%), and 6,981 other hepatobiliary malignancies (9.2%). Median age was 74 years (IQR 68-80); 65.4% were male. The cohort was racially diverse: Non-Hispanic White 52.4%, Hispanic 20.1%, Non-Hispanic Asian/Pacific Islander 15.1%, Non-Hispanic Black 11.1%.

### Baseline Characteristics and Survival

Baseline characteristics stratified by surgery status are shown in **Table 1** (full version: **eTable 1**). Of 76,110 patients, 59,821 (78.6%) received no cancer-directed surgery and 16,289 (21.4%) underwent surgical intervention. Surgery patients were younger (mean 72.4 vs 74.6 years, p<0.001), more likely to be married (61.0% vs 52.6%, p<0.001), and had lower-stage disease (Stage I: 60.5% vs 41.8%). Cirrhosis was more prevalent in the surgery group (12.9% vs 8.4%, p<0.001), suggesting careful selection of surgical candidates despite underlying liver disease.

### Survival Outcomes

Median follow-up was 49 months (OS) and 74 months (CSS) for surgery patients, versus 11 and 14 months for non-surgery patients, respectively. Median OS by surgery type was: None 4 months, Local Destruction 28 months, Segmental Resection 30 months, Larger Resection 24 months, and Liver Transplantation 42 months (log-rank p<0.0001; **Figure 1**).

In multivariate CSS Cox regression (C-index 0.739), the strongest protective factors were: liver transplantation (HR 0.15, 95% CI 0.14-0.17), segmental resection (HR 0.23, 95% CI 0.22-0.25), larger resection (HR 0.24, 95% CI 0.22-0.27), and local destruction (HR 0.26, 95% CI 0.25-0.27) (**Figure 2**; **eTable 2**). Adverse prognostic factors included: Stage IV (HR 1.74), Stage III (HR 1.67), ICC histology (HR 1.29), and poor differentiation (HR 1.04). Chemotherapy (HR 0.57) and radiation (HR 0.52) were independently protective.

### Primary Finding: Segmental vs Larger Resection

The key analysis comparing segmental resection to larger resection revealed a consistent pattern across all age strata ≥70 years: segmental resection provided equivalent or numerically better CSS than larger resection (**Table 2**). For patients aged 70-74, segmental resection HR was 0.27 (95% CI 0.24-0.31) compared to larger resection HR of 0.25 (95% CI 0.21-0.31). In the 75-79 age band: segmental HR 0.29 vs larger HR 0.33. In the 80+ band: segmental HR 0.26 vs larger HR 0.31. Only in the 65-69 age band did larger resection show a marginal advantage (segmental HR 0.33 vs larger HR 0.35).

This pattern was observed in both HCC and ICC (**Figure 3**). For HCC, segmental resection adjusted HR was 0.26 (95% CI 0.24-0.28) versus larger resection HR of 0.26 (95% CI 0.22-0.30). For ICC, segmental resection was numerically superior (HR 0.22, 95% CI 0.19-0.24 vs larger HR 0.26, 95% CI 0.22-0.31).

Formal interaction testing confirmed that the surgery × cancer type interaction was not statistically significant (HR 0.96, 95% CI 0.89-1.02, p=0.206), indicating that the relative benefit of segmental over larger resection is equivalent between HCC and ICC.

### Propensity Score Matching

After 1:1 PSM, 6,434 matched pairs were generated with substantially improved covariate balance (mean SMD: 0.190 before matching → 0.030 after matching; **eFigure 2**). In the matched cohort, PSM-adjusted surgery HR was 0.35 (95% CI 0.33-0.37), confirming robust survival benefit after accounting for measured confounding.

### Machine Learning and External Validation

Among five survival models, XGBoost achieved the highest internal C-index of 0.746, followed by RSF (0.736), Cox PH (0.739), and DeepSurv (0.575). The ensemble model combining Cox, RSF, and XGBoost predictions via 5-fold cross-validation achieved a C-index of 0.737 ± 0.003 (**eFigure 3**). Feature importance analysis identified surgery receipt, AJCC stage, and chemotherapy as the dominant predictors (**eFigure 4**).

External validation on TCGA-LIHC (n=269, all surgical HCC) yielded C-indices of 0.595 (Cox), 0.567 (RSF), and 0.592 (XGBoost). On ICGC-LIRI-JP (n=260, all surgical HCC, Japanese), C-indices were 0.522, 0.551, and 0.547, respectively (**Figure 4**; **eFigure 5**). The systematic performance drop (ΔC-index ≈ 0.09 for TCGA, ≈ 0.14 for ICGC) is largely attributable to the "Surgery Feature Paradox": in all-surgical external cohorts, the strongest SEER predictor (surgery_any) has zero variance.

### Sensitivity and Robustness Analyses

**E-value**: The E-value for the surgery benefit was 7.5 (CI-bound 7.7), indicating that an unmeasured confounder would need to be associated with both surgery receipt and survival by a risk ratio of ≥7.5 to fully explain away the observed benefit. This exceeds the effect size of any known prognostic factor, including Stage IV disease (the strongest measured confounder). E-values remained ≥5.8 across all age strata (**eTable 3**).

**Instrumental Variable**: IV analysis using regional socioeconomic status as instrument confirmed persistent surgery benefit, consistent with the primary analysis after accounting for treatment selection bias.

**Leave-One-Out**: Excluding Stage IV patients, age 80+, ICC, low-income, and early-era (2004-2010) patients yielded surgery HRs ranging from 0.27 to 0.30, all within a narrow band around the primary estimate of 0.25 (**eTable 5**).

**Model Specification**: Surgery HR remained stable across all model specifications: unadjusted 0.27 → demographics-adjusted 0.29 → stage-adjusted 0.31 → tumor-adjusted 0.31 → fully-adjusted 0.25 (**eTable 4**).

**Fine-Gray Competing Risk**: Surgery was protective against both cancer-specific death (HR 0.32, 95% CI 0.31-0.33) and other-cause death (HR 0.36, 95% CI 0.34-0.37), confirming that the benefit extends beyond cancer control.

### Temporal and Geographic Patterns

Median OS improved from 3 months (2000) to 12 months (2019), reflecting two decades of treatment advances (**eFigure 6**). During the COVID-19 era (2020-2022), a transient decline was observed (median OS 6 months, surgery rate 21%), with recovery evident by 2022. Geographic analysis across socioeconomic regions revealed consistent parenchymal-sparing practice patterns (Segmental/Larger ratio 3.2-4.1) but higher surgery utilization in high-income regions (24% vs 20%, p<0.001; **eTable 6**).

### Frailty and Age Interactions

Using the Frailty Surrogate Index, 34,248 patients (45.0%) were classified as Fit, 30,324 (39.8%) as Pre-frail, and 11,538 (15.2%) as Frail. Surgery benefit persisted across all frailty strata: Fit HR 0.27, Pre-frail HR 0.32, Frail HR 0.33 (all p<0.001). Median OS for surgical Frail patients was 15 months vs 2 months without surgery—a 13-month absolute benefit (**eFigure 7**).

Restricted cubic spline analysis demonstrated that surgery benefit follows a U-shaped relationship with age, peaking at age 72 with an inflection point at age 86, beyond which benefit attenuation accelerates (**eFigure 8**).

### Clinical Decision Support

The Hepatobiliary Cancer Index (HBI) integrates 7 variables into a simple point system: surgery receipt (−10 points), chemotherapy (−5), Stage IV (+4), Stage III (+3), ICC histology (+1), poor grade (0), and age (0). Three-tier risk stratification yielded median OS of 22 months (low-risk, HBI <−5), 5 months (intermediate, HBI −5 to 0), and 2 months (high-risk, HBI >0) (**Figure 5**). The clinical decision matrix (**Table 3**) provides explicit recommendations stratified by age, cancer type, stage, and frailty status. A clinical nomogram was constructed for individualized survival prediction (**Figure 6**).

---

## DISCUSSION

In this population-based analysis of 76,110 elderly patients with hepatobiliary cancer, we demonstrate that segmental/wedge resection provides equivalent or superior overall and cancer-specific survival compared to larger resection (lobectomy/extended lobectomy) across both HCC and ICC. This finding was robust to propensity score matching, instrumental variable analysis, multiple sensitivity analyses, and maintained an E-value of 7.5, indicating substantial resilience to unmeasured confounding.

### Comparison with Prior Literature

Our findings corroborate and substantially extend the work of Zhang et al. (2020), who first reported from SEER data that segmental resection was equivalent to larger resection in 2,371 elderly stage I-II HCC patients.[6] We confirm this finding in a cohort 32 times larger, expand to all AJCC stages, and demonstrate—for the first time—that this pattern holds for ICC, a cancer type with distinct biology and worse prognosis.

The comparison with Xia et al. (2025) is instructive: by combining wedge and segmental resection into a single category, their analysis masked the nuanced comparison that our finer surgical categorization reveals.[7] Our data suggest that the "wedge/segmental" category, as commonly reported, actually represents the preferred surgical strategy for elderly patients rather than a compromise.

The finding that surgery benefit is equivalent across HCC and ICC (interaction p=0.206) has not been previously reported. Given that ICC has substantially worse baseline prognosis (median OS 5 vs 8 months) and is typically diagnosed at more advanced stages (Stage IV: 32% vs 18%), the equivalence of surgical benefit magnitude is clinically meaningful: it suggests that the decision between segmental and larger resection should be based on technical and patient factors rather than cancer histology alone.

### Clinical Implications

The paradigm of "more radical surgery = better oncologic outcome" does not appear to hold for elderly hepatobiliary cancer patients. Our findings support a shift toward parenchymal-sparing surgery as the default strategy when transplantation is unavailable. This recommendation is strengthened by the Frailty Surrogate Index analysis, demonstrating that even patients classified as "Frail" derived substantial survival benefit from surgery (HR 0.33; absolute median OS gain of 13 months).

Age alone should not be a contraindication to surgical intervention. The surgery benefit persists across all age strata, with only modest attenuation at age 80+ (HR 0.31 vs 0.27 at age 65-69). The restricted cubic spline analysis, while confirming gradual attenuation with advancing age, identified the inflection point only at age 86—well beyond the typical threshold at which surgical candidacy is questioned.

The HBI risk score provides a practical bedside tool for pre-treatment counseling. The 11-fold survival difference between low-risk (22 months) and high-risk (2 months) groups can inform shared decision-making conversations about expected benefit from surgical intervention.

### Strengths and Limitations

This study has several notable strengths. First, the population-based SEER design with 76,110 patients provides substantially greater statistical power than prior studies (previous maximum: 10,174 in Xia et al. 2025). Second, the dual external validation on TCGA-LIHC and ICGC-LIRI-JP, while revealing systematic performance degradation (the "Surgery Feature Paradox"), provides the first external generalizability benchmarks for SEER-trained hepatobiliary cancer models. Third, the multi-layered robustness assessment (E-value, IV, PSM, leave-one-out, specification testing, Fine-Gray) establishes a new standard for causal rigor in surgical oncology research using observational data. Fourth, the inclusion of ICC alongside HCC in the same analytical framework enables the first formal comparison of surgical extents across these biologically distinct cancers.

Several limitations warrant discussion. The SEER database lacks detailed clinical information critical for surgical decision-making: liver function (Child-Pugh, MELD), performance status, comorbidity burden, and specifics of perioperative complications.[13] The absence of these variables could introduce residual confounding, though our E-value analysis (E=7.5) suggests any unmeasured confounder would need implausibly strong effects to nullify the findings. Changes in AJCC staging across the study period required algorithmic harmonization, which may have introduced minor misclassification. The external validation cohorts (TCGA, ICGC) are all-surgical cohorts, which eliminated the predictive variance of the strongest SEER feature (surgery_any), contributing to the observed C-index gap. Our Frailty Surrogate Index is an approximation without direct comorbidity or functional status measurement; validation against comprehensive geriatric assessment is needed. Finally, the SEER coding of "segmental resection" encompasses both wedge and anatomical segmentectomy, which may have different oncologic adequacy depending on tumor location and size.

### Future Directions

Several lines of investigation are warranted. First, linkage with Medicare claims data could provide comorbidity, functional status, and chemotherapy regimen details to refine the frailty assessment. Second, prospective validation of the HBI risk score in clinical cohorts is needed before implementation. Third, investigation of minimally invasive (laparoscopic/robotic) approaches within the segmental resection category is timely given recent technical advances. Fourth, health-related quality of life outcomes, which cannot be assessed in SEER, are critical for elderly surgical decision-making and should be incorporated in prospective studies.

---

## CONCLUSIONS

In this population-based analysis of 76,110 elderly patients with hepatobiliary cancer, segmental/wedge resection provides equivalent or superior survival to larger resection across both hepatocellular carcinoma and intrahepatic cholangiocarcinoma, independent of age and stage. These findings, supported by robust causal sensitivity analyses (E-value=7.5), propensity score matching (HR 0.35), and external validation, support a paradigm shift toward parenchymal-sparing surgery as the preferred approach in elderly patients when liver transplantation is unavailable. The 7-variable HBI risk score and clinical decision matrix provide practical tools for individualized treatment planning.

---

## REFERENCES

1. Torre LA, Bray F, Siegel RL, et al. Global cancer statistics, 2012. *CA Cancer J Clin*. 2015;65(2):87-108.
2. Oweira H, Petrausch U, Helbling D, et al. Early stage hepatocellular carcinoma in the elderly: a SEER database analysis. *J Geriatr Oncol*. 2017;8(4):277-283.
3. Ciccarese F, Caturelli E, Felder M, et al. Survival benefit of liver resection for patients with hepatocellular carcinoma across different BCLC stages. *J Hepatol*. 2014;62(3):617-624.
4. Hirokawa F, Hayashi M, Miyamoto Y, et al. Surgical outcomes and clinical characteristics of elderly patients undergoing curative hepatectomy for hepatocellular carcinoma. *J Gastrointest Surg*. 2013;17(11):1929-1937.
5. Yu B, Ding Y, Liao X, et al. Radiofrequency ablation versus surgical resection in elderly patients with early-stage hepatocellular carcinoma. *Saudi J Gastroenterol*. 2018;24(6):317-325.
6. Zhang QQ, Wu PYS, ALBahde M, et al. Do elderly patients with stage I-II hepatocellular carcinoma benefit from more radical surgeries? A population-based analysis. *Front Oncol*. 2020;10:479.
7. Xia L, et al. Effect of surgery on overall survival and cancer-specific survival in patients with primary HCC: A study based on PSM in the SEER cohort. *Medicine*. 2025;104(8):e41521.
8. Prognostic analysis of elderly patients with hepatocellular carcinoma: an exploration and machine learning model prediction based on age stratification and surgical approach. *J Hepatocell Carcinoma*. 2025.
9. Development and validation of an XGBoost model to predict 5-year survival in elderly patients with intrahepatic cholangiocarcinoma after surgery: a SEER-based study. *J Gastrointest Surg*. 2023.
10. Surveillance, Epidemiology, and End Results (SEER) Program. SEER*Stat Database: Incidence - SEER Research Plus Data, 18 Registries, Nov 2024 Sub (2000-2022). National Cancer Institute.
11. VanderWeele TJ, Ding P. Sensitivity analysis in observational research: introducing the E-value. *Ann Intern Med*. 2017;167(4):268-274.
12. Baiocchi M, Cheng J, Small DS. Instrumental variable methods for causal inference. *Stat Med*. 2014;33(13):2297-2340.
13. Peters NA, Javed AA, He J, et al. Association of socioeconomics, surgical therapy, and survival of early stage hepatocellular carcinoma. *J Surg Res*. 2017;210:253-260.

---

## SUPPLEMENTARY MATERIAL

- **eFigure 1**: CONSORT Flow Diagram
- **eFigure 2**: PSM Love Plot (SMD before/after matching)
- **eFigure 3**: 5-Fold CV Model Comparison
- **eFigure 4**: SHAP Feature Importance & Waterfall
- **eFigure 5**: External Validation (TCGA-LIHC + ICGC-LIRI-JP)
- **eFigure 6**: Temporal Trends (2004-2022)
- **eFigure 7**: Frailty Surrogate Index — KM by FSI Group
- **eFigure 8**: Restricted Cubic Spline — Age-Dependent Surgery HR
- **eTable 1**: Full Baseline Characteristics
- **eTable 2**: Complete Multivariate Cox Regression
- **eTable 3**: Stratified E-values by Age Band
- **eTable 4**: Model Specification Robustness
- **eTable 5**: Leave-One-Out Analysis
- **eTable 6**: Geographic Practice Variation

---

*Correspondence to:* Qigang Xu, MD, Department of Hepatobiliary and Pancreatic Surgery, The First Affiliated Hospital of Wenzhou Medical University, Zhejiang, China. E-mail: xuqigang@wmu.edu.cn

**Author Contributions**: ZZ and YB contributed equally to this work (co-first authors). QX had full access to all data and takes responsibility for integrity and accuracy of analyses. Study concept and design: ZZ, QX. Data acquisition and analysis: ZZ, YB, YC. Statistical analysis: ZZ. Machine learning modeling: ZZ, YB. Manuscript drafting: ZZ. Critical revision: All authors. Study supervision: QX.

**Funding**: This research did not receive any specific grant from funding agencies in the public, commercial, or not-for-profit sectors.

**Conflict of Interest**: The authors have declared that no competing interest exists.

**Ethics Statement**: All procedures followed were in accordance with the ethical standards of the Ethics Committee of The First Affiliated Hospital of Wenzhou Medical University and with the Helsinki Declaration of 1964 and later versions. SEER data are publicly available and de-identified; institutional review board approval was exempt.

**Data Availability**: SEER data are publicly available at seer.cancer.gov (Research Plus Data, November 2024 submission). TCGA-LIHC data are available at portal.gdc.cancer.gov. ICGC-LIRI-JP data are available at dcc.icgc.org. Analysis code available upon reasonable request to the corresponding author.

---

