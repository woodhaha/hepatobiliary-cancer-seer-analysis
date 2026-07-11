# Submission Package — Hepatobiliary Cancer Manuscript

## Target Journal: Annals of Surgical Oncology (Primary) / JAMA Surgery (Stretch)

---

## Required Materials Checklist

| Item | Status | File |
|------|--------|------|
| Title Page | ✅ | Included in manuscript |
| Manuscript Body | ✅ | `05_Submission/manuscript.md` |
| Abstract (structured) | ✅ | 250 words |
| Key Words | ✅ | 6 keywords |
| Figures (main) | ✅ | Fig1-6 (6 main + 8 suppl) |
| Tables (main) | ✅ | Table 1-3 |
| Supplementary Materials | ✅ | `supplementary_materials.md` |
| Cover Letter | ⚠️ | Draft below |
| Author Contributions | ✅ | Included |
| COI Statement | ✅ | No competing interests |
| Funding Statement | ✅ | No specific funding |
| Ethics Statement | ✅ | Included |
| Data Availability | ✅ | SEER public + TCGA/ICGC public |
| Title Page with Author Info | ✅ | Included in manuscript header |

---

## DRAFT Cover Letter

```
Date: [Insert Date]

Dear Editor,

We are pleased to submit our manuscript entitled "Parenchymal-Sparing vs Major
Hepatectomy in Elderly Hepatobiliary Cancer: A Population-Based Analysis with
Machine Learning Ensemble and External Validation" for consideration in Annals of
Surgical Oncology.

This study addresses a critical clinical question: whether parenchymal-sparing
liver resection (segmental/wedge) provides equivalent oncologic outcomes to major
hepatectomy in elderly patients with hepatobiliary cancer. Using SEER data covering
76,110 patients, we demonstrate for the first time that segmental resection is
equivalent or superior to larger resection across both hepatocellular carcinoma and
intrahepatic cholangiocarcinoma—a finding with direct implications for surgical
decision-making in an aging population.

Key strengths of this manuscript include:
1. Population-based design with 76,110 patients—the largest analysis of its kind
2. First simultaneous comparison of surgical extents in HCC and ICC
3. Robust causal sensitivity analysis (E-value=5.2, PSM, instrumental variable)
4. Machine learning ensemble with 5-fold cross-validation
5. Dual external validation (TCGA-LIHC + ICGC-LIRI-JP)
6. Practical clinical decision matrix and risk stratification tool

All authors have read and approved the final manuscript. The study used publicly
available de-identified data and was exempt from institutional review board approval.
The authors declare no competing interests.

We believe this manuscript aligns well with Annals of Surgical Oncology's focus on
clinical outcomes research in surgical oncology and will be of significant interest
to your readership.

Sincerely,

Zhuha Zhou, MD
Corresponding Author
Department of Gastroenterology Surgery
The First Affiliated Hospital of Wenzhou Medical University
Zhejiang, China
E-mail: zhouzhuha@wmu.edu.cn
```

---

## Figure List

### Main Figures

| Figure | Title | File |
|--------|-------|------|
| **Fig 1** | Kaplan-Meier Survival Curves (6 panels: OS/CSS by surgery, stage, age, cancer type) | `Fig1_KM.pdf` |
| **Fig 2** | Multivariate Cox Forest Plot (18 variables) | `Fig2_Forest.pdf` |
| **Fig 3** | HCC vs ICC Surgical Strategy Comparison | `Fig20_HCCvsICC.pdf` |
| **Fig 4** | External Validation (SEER→TCGA+ICGC) | `Fig5_ExternalValidation.pdf` |
| **Fig 5** | Composite Analysis (CV + PSM + HBI + Age-Surgery) | `Fig6_CompositeAnalysis.pdf` |
| **Fig 6** | Clinical Nomogram | `Fig10_Nomogram.png` |

### Supplementary Figures (eFigure 1-8)

| eFigure | Title | File |
|---------|-------|------|
| **eFig 1** | CONSORT Flow Diagram | `Fig19_CONSORT.pdf` |
| **eFig 2** | PSM Love Plot | `Fig6_CompositeAnalysis.pdf` (Panel B) |
| **eFig 3** | 5-Fold CV Model Comparison | `Fig6_CompositeAnalysis.pdf` (Panel A) |
| **eFig 4** | SHAP Feature Importance | `FigS1_SHAP.png` |
| **eFig 5** | External Validation Details | `Fig5_ExternalValidation.pdf` |
| **eFig 6** | Temporal Trends (2004-2022) + COVID | `FigS4_TemporalTrends.png` + `Fig15_COVID.png` |
| **eFig 7** | Frailty Surrogate Index | `Fig22_Frailty.png` |
| **eFig 8** | Age-Surgery RCS + Landmark + Competing Risk | `Fig21_RCS_AgeSpline.png` + `FigS2_CompetingRisk.pdf` |

### Additional Supplementary

| Item | File |
|------|------|
| Supplementary eTable S1-S6 | `supplementary_materials.md` |
| Time-Dependent AUC | — *(not generated)* |
| LASSO Variable Selection | `Fig7_LASSO.png` |
| Subgroup Forest Plot | `Fig11_Subgroup.png` |
| Landmark Conditional Survival | `Fig12_Landmark.png` |
| AFP Analysis | `Fig13_AFP.png` |
| Race Disparities | `Fig14_Race.png` |
| SES Interaction | — *(not generated)* |
| Geographic Variation | `Fig23_Geography.png` |

---

## Table List

### Main Tables

| Table | Title | Location |
|-------|-------|----------|
| **Table 1** | Baseline Characteristics | Manuscript body + eTable 1 |
| **Table 2** | Age-Stratified Surgery CSS Hazard Ratios (HCC vs ICC) | Manuscript body |
| **Table 3** | Clinical Decision Matrix | Manuscript body |

### Supplementary Tables

| eTable | Title |
|--------|-------|
| eTable 1 | Full Baseline Characteristics |
| eTable 2 | Complete Multivariate Cox Regression |
| eTable 3 | Stratified E-values by Age Band |
| eTable 4 | Model Specification Robustness |
| eTable 5 | Leave-One-Out Sensitivity Analysis |
| eTable 6 | Geographic Practice Variation |

---

## Suggested Reviewers

1. **Wei-Lin Wang** (Zhejiang University) — Published the reference SEER-HCC paper (Zhang et al. 2020, Front Oncol)
2. **San-Lin You** (Fu Jen Catholic University) — Reviewed the reference paper; expertise in HCC epidemiology
3. **Lei Xia**, MD — Corresponding author of Xia et al. 2025 (*Medicine*), the most recent SEER-based HCC surgical comparison
4. **Expert in surgical oncology with SEER methodology**
5. **Biostatistician with survival ML expertise**

---

## Opposed Reviewers

None.
