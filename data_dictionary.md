# Data Dictionary — Hepatobiliary Cancer SEER Analysis

> Generated: 2026-06-26 · Source: SEER Research Plus Data (Nov 2024)

## File: 02_Data/raw/hepatobiliary cancer.csv

| # | Variable | Type | Description | Coding |
|---|----------|------|-------------|--------|
| 1 | Patient ID | int | SEER patient identifier | — |
| 2 | Type of Reporting Source | cat | Data source type | Hospital inpatient/outpatient, Laboratory only, etc. |
| 3 | Marital status at diagnosis | cat | Marital status | Married, Divorced, Single, Widowed, etc. |
| 4 | Median household income | cat | Census tract income, inflation adj to 2023 | Binned categories |
| 5 | Rural-Urban Continuum Code | cat | RUCC urbanization level | 1-9 |
| 6 | Sex | cat | Biological sex | Male, Female |
| 7 | Year of diagnosis | int | Calendar year of diagnosis | 2000-2022 |
| 8 | PRCDA 2020 | cat | PRCDA region indicator | — |
| 9 | Race recode | cat | Race (W, B, AI, API) | White, Black, American Indian, Asian/Pacific Islander |
| 10 | Origin recode NHIA | cat | Hispanic origin | Hispanic, Non-Hispanic |
| 11 | Race and origin recode | cat | Combined race+ethnicity | NHW, NHB, NHAIAN, NHAPI, Hispanic |
| 12 | Age recode with <1 and 90+ | cat | Age group binned | 45-49, 50-54, ..., 85+ |
| 13 | Site recode ICD-O-3 | cat | Primary site | C22.0-Liver, C22.1-Intrahepatic Bile Duct, etc. |
| 14 | Histologic Type ICD-O-3 | cat/int | Histology code | 8170 = HCC, 8160 = Cholangiocarcinoma |
| 15 | Grade Recode (thru 2017) | cat | Differentiation grade | Well, Moderately, Poorly, Undifferentiated, Unknown |
| 16 | Derived Summary Grade 2018 | cat | Grade (2018+) | 1-9 |
| 17 | Laterality | cat | Tumor laterality | Not a paired site, etc. |
| 18 | Diagnostic Confirmation | cat | Diagnostic method | Positive histology, Radiography, etc. |
| 19 | Combined Summary Stage | cat | SEER summary stage | Localized, Regional, Distant, etc. |
| 20 | AJCC Stage 6th ed (2004-2015) | cat | AJCC TNM stage 6th | I, II, III, IV, UNK |
| 21-24 | Derived AJCC T/N/M 6th | cat | T, N, M stage 6th ed | TX-T4, NX-N1, MX-M1 |
| 25-28 | Derived AJCC Stage/T/N/M 7th | cat | Stage 7th ed (2010-2015) | — |
| 29-32 | Derived EOD 2018 Stage/T/N/M | cat | EOD staging (2018+) | — |
| 33 | RX Summ--Surg Prim Site | cat | Primary site surgery | 0=None, 30=Segmental/Wedge, 40=Lobectomy, 50=Extended lobectomy, 60=Transplant |
| 34 | RX Summ--Scope Reg LN Sur | cat | Lymph node surgery scope | 0=None, 1-4=Regional LN |
| 35 | Reason no cancer-directed surgery | cat | Reason for no surgery | Not recommended, Recommended but refused, etc. |
| 36 | Radiation recode | cat | Radiation therapy | None, Beam radiation, etc. |
| 37 | Chemotherapy recode | cat | Chemotherapy | Yes, No/Unknown |
| 38 | Time from diagnosis to treatment | cat | Days to treatment | Binned categories |
| 39-42 | SEER Combined Mets at DX | bin | Metastasis sites (bone/brain/liver/lung) | No, Yes |
| 43 | COD to site recode | cat | Cause of death site | Liver, Intrahepatic Bile Duct, etc. |
| 44 | SEER cause-specific death | cat | CSS death classification | Dead (attributable), Alive or dead of other |
| 45 | SEER other cause of death | cat | Other cause death classification | — |
| 46 | Survival months | float | Follow-up months | 0-XXX |
| 47 | Vital status recode | cat | Vital status at cutoff | Alive, Dead |
| 48 | Age recode with single ages | int | Exact age at diagnosis | 0-99 |
| 49 | CS tumor size (2004-2015) | int | Tumor size in mm | 1-989 (999=Unknown) |
| 50 | Regional nodes examined | int | # nodes examined | 0-XX (98=Unknown) |
| 51 | Regional nodes positive | int | # nodes positive | 0-XX (98=Unknown) |
| 52 | Primary Site - labeled | cat | Primary site label | C22.0-Liver, etc. |
| 53 | Total # in situ/malignant tumors | int | Total tumor count | 1-XX |
| 54 | AFP Pretreatment Interpretation | cat | Alpha-fetoprotein | Negative, Borderline, Positive, Unknown |
| 55 | Fibrosis Score | cat | Ishak fibrosis score | 0-4, 5-6, Unknown |

## Derived Variables (to be engineered)

| Variable | Type | Derivation | Source |
|----------|------|------------|--------|
| age_c | float | Age centered at 67 | Age recode with single ages |
| elderly | bin | Age >= 65 | age_c >= -2 |
| surgery_any | bin | Any surgical resection | RX Summ--Surg Prim Site > 0 |
| surgery_type | cat | None/Local/Segmental/Larger/Transplant | Recode of surgery codes |
| stage_I/II/III/IV | bin | AJCC stage dummies | Combined staging |
| hcc | bin | HCC vs ICC | Site + Histology |
| os_status | bin | All-cause death | Vital status recode |
| css_status | bin | Cancer-specific death | SEER cause-specific death |
| os_time | int | Overall survival months | Survival months |
| css_time | int | CSS follow-up months | Survival months (censored for other death) |

## SEER Surgery Codes Mapping

| Code | Category | Definition |
|------|----------|------------|
| 0 | None | No surgery |
| 10-27 | Local Destruction | RFA, PEI, cryoablation |
| 30 | Segmental Resection | Wedge/segmental resection |
| 40 | Larger Resection | Lobectomy |
| 50 | Larger Resection | Extended lobectomy |
| 60 | Liver Transplantation | Total hepatectomy + transplant |
