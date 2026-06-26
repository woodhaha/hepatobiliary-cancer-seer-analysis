# Clinical Information Comparison — All Data Sources

## Cohort Overview

| Variable | SEER | TCGA-LIHC | ICGC-LIRI-JP |
|----------|------|-----------|-------------|
| **Source** | SEER 18 Registries | TCGA (US multi-center) | ICGC (Japan) |
| **Years** | 2004–2022 | 2013–2017 | 2005–2014 |
| **Total N** | 76,110 (elderly) | 269 | 260 |
| **HCC** | 57,380 (75.4%) | 269 (100%) | 260 (100%) |
| **ICC** | 11,749 (15.4%) | 0 | 0 |
| **Other** | 6,981 (9.2%) | 0 | 0 |
| **Population** | US (28% coverage) | US multi-center | Japan |
| **Data type** | Population registry | Genomic + clinical | Genomic + clinical |

---

## Demographics

| Variable | SEER | TCGA | ICGC |
|----------|------|------|------|
| **Median Age** | 74 (IQR 68–80) | 59 (range 19–90) | 67 (range 31–89) |
| **Age ≥65** | 100% (by design) | ~35% | ~50% |
| **Male %** | 65.4% | 67% | 72% |
| **Race** | NHW 52%, Hisp 20%, NHAPI 15%, NHB 11% | NHW 52%, Asian 28%, NHB 5% | Asian 100% |
| **Married %** | 54.4% | — | — |
| **Income data** | Census tract median | — | — |

---

## Tumor Characteristics

| Variable | SEER | TCGA | ICGC |
|----------|------|------|------|
| **AJCC Stage I** | 45.8% | 52% (pathologic) | 24% |
| **Stage II** | 16.8% | 29% | 37% |
| **Stage III** | 14.9% | 19% | 30% |
| **Stage IV** | 22.5% | 0% (all resectable) | 9% (resectable) |
| **Grade (Poor)** | 31.6% | ~25% (G3/G4) | — |
| **Tumor Size** | Median 75mm (SEER: CS Tumor Size) | — | — |
| **Histology** | HCC 75%, ICC 15% | HCC 100% | HCC 100% |
| **AFP available** | Yes (43% missing) | Yes (serum) | — |
| **Cirrhosis** | 9.4% (documented; clinical reality ~80%) | 25% (Ishak fibrosis score) | — |

---

## Treatment

| Variable | SEER | TCGA | ICGC |
|----------|------|------|------|
| **Surgery %** | 21.4% | 100% (all resection) | 100% (all resection) |
| **Segmental** | 3.3% | Unknown | Unknown |
| **Larger** | 0.9% | Unknown | Unknown |
| **Transplant** | 2.3% | 0% | 0% |
| **Local Destruction** | 14.2% | 0% | 0% |
| **Chemotherapy** | 31.9% | ~10% (neoadjuvant) | — |
| **Radiation** | 12.7% | ~5% | — |
| **Laparoscopic vs Open** | Not available | Not available | Not available |

---

## Liver-Specific Variables

| Variable | SEER | TCGA | ICGC |
|----------|------|------|------|
| **Cirrhosis / Fibrosis** | Ishak score (51% missing) | Ishak score (path review) | — |
| **Child-Pugh** | **NOT available** | Available | — |
| **MELD** | **NOT available** | — | — |
| **Portal hypertension** | **NOT available** | — | — |
| **BCLC stage** | **NOT available** (only AJCC TNM) | — | — |
| **Viral hepatitis** | Not in SEER | HBV/HCV serology | HBV-dominant |
| **AFP** | Available (2010+) | Available | — |
| **Albumin** | **NOT available** | Available | — |
| **Bilirubin** | **NOT available** | Available | — |
| **Vascular invasion** | Not in SEER | Available (path) | — |

---

## Survival & Outcomes

| Variable | SEER | TCGA | ICGC |
|----------|------|------|------|
| **Survival measure** | OS + CSS | OS | OS |
| **Median OS** | 8m (HCC), 5m (ICC) | 22m | 21m |
| **Dead %** | 82.7% | 7% | 32% |
| **Cancer-specific death** | 68.5% | — | — |
| **Median follow-up** | 49m (OS), 74m (CSS) | ~24m | ~24m |
| **Max follow-up** | 275 months | ~120 months | ~120 months |
| **Censoring rate** | 17.3% | 93% | 68% |

---

## Variables NOT Available in Any Source

| Variable | SEER | TCGA | ICGC | Importance |
|----------|------|------|------|------------|
| **Performance Status (ECOG)** | ✗ | ✗ | ✗ | Critical for elderly |
| **Comorbidity (CCI/ACCI)** | ✗ | ✗ | ✗ | Critical for elderly |
| **Frailty (clinical)** | ✗ | ✗ | ✗ | Our FSI is proxy only |
| **Resection margin (R0/R1/R2)** | ✗ | ✗ | ✗ | Major for ICC |
| **Lymph node dissection** | Partial | ✗ | ✗ | Quality metric for ICC |
| **Perioperative complications** | ✗ | ✗ | ✗ | Clavien-Dindo grading |
| **90-day mortality** | ✗ | ✗ | ✗ | Surgical quality |
| **Sarcopenia (CT-SMI)** | ✗ | ✗ | ✗ | Novel prognostic factor |
| **NLR/PLR (inflammation)** | ✗ | ✗ | ✗ | Novel prognostic factor |

---

## External Validation Feasibility Map

| Variable group | SEER → TCGA? | SEER → ICGC? | Why gap? |
|---------------|-------------|-------------|----------|
| Age | Partial | Partial | TCGA younger, ICGC older |
| Sex | ✓ | ✓ | Consistent |
| Race | Partial | ✗ | ICGC all Asian |
| Stage | ✓ | Partial | Different staging calibration |
| Surgery type | ✗ | ✗ | All-surgical = zero variance |
| Chemotherapy | Partial | — | Different indications |
| Survival | Partial | Partial | Very different baseline risk |
| Liver function | ✗ | ✗ | Not in SEER at all |
| **Overall** | **0.595** | **0.547** | Surgery paradox + population drift |

---

## Key Takeaway

SEER provides the largest sample but lacks the clinical variables (Child-Pugh, performance status) that drive real-world surgical decision-making. TCGA and ICGC provide richer tumor biology (genomics, pathology) but are all-surgical, younger, and HCC-only — limiting their utility as external validators for an elderly mixed-histology treatment recommendation model. **The missing variables (ECOG, CCI, Child-Pugh) affect all three sources equally**, making this a universal limitation in population-based hepatobiliary surgical research.
