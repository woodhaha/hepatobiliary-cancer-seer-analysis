# Hepatobiliary Cancer Survival Analysis

**Cohort**: 76110 elderly (≥65), 59821 Non-surgery, 16289 Surgery

## Table 1: Baseline Characteristics

| Variable | Non-surgery (N=59821) | Surgery (N=16289) | P-value |
|----------|-------------------|----------------|---------|
| Age (mean±SD) | 74.6±7.0 | 72.4±5.8 | <0.001 |
| Male (%) | 65.2% | 65.8% | 0.155 |
| Married (%) | 52.6% | 61.0% | <0.001 |
| Stage I (%) | 41.8% | 60.5% | 0.500 |
| Stage II (%) | 16.1% | 19.6% | 0.500 |
| Stage III (%) | 16.1% | 10.3% | 0.500 |
| Stage IV (%) | 26.0% | 9.7% | 0.500 |
| Poor Grade (%) | 32.7% | 27.4% | <0.001 |
| ICC (%) | 16.0% | 13.5% | <0.001 |
| Chemotherapy (%) | 33.3% | 26.9% | <0.001 |
| Radiation (%) | 14.4% | 6.6% | <0.001 |
| Cirrhosis (%) | 8.4% | 12.9% | <0.001 |

## Figure 1: Kaplan-Meier

✓ Fig1 saved

## Log-Rank Tests

OS Surgery vs Non: χ²=13050.7, p=0.0000
CSS Surgery vs Non: χ²=11716.9, p=0.0000

## Cox Regression

### Univariate CSS Cox
| Variable | HR | 95% CI | P |
|---|---|---|---|
| age_c | 1.03 | [1.03-1.03] | 0.000 |
| male | 0.96 | [0.94-0.98] | 0.000 |
| married | 0.87 | [0.86-0.89] | 0.000 |
| race_nhb | 1.01 | [0.98-1.04] | 0.522 |
| race_nhapi | 0.77 | [0.75-0.79] | 0.000 |
| race_hispanic | 1.00 | [0.97-1.02] | 0.671 |
| stage_2 | 0.85 | [0.83-0.87] | 0.000 |
| stage_3 | 1.45 | [1.41-1.48] | 0.000 |
| stage_4 | 2.03 | [1.99-2.07] | 0.000 |
| grade_poor | 1.05 | [1.03-1.07] | 0.000 |
| is_icc | 1.36 | [1.33-1.39] | 0.000 |
| surg_local_destruction | 0.32 | [0.31-0.33] | 0.000 |
| surg_segmental_resection | 0.39 | [0.37-0.42] | 0.000 |
| surg_larger_resection | 0.42 | [0.37-0.46] | 0.000 |
| surg_transplant | 0.20 | [0.19-0.22] | 0.000 |
| chemotherapy | 0.72 | [0.71-0.73] | 0.000 |
| radiation | 0.68 | [0.66-0.70] | 0.000 |
| cirrhosis | 0.63 | [0.61-0.65] | 0.000 |

### Multivariate CSS Cox
| Variable | HR | 95% CI | P |
|---|---|---|---|
| age_c | 1.02 | [1.01-1.02] | <0.001 |
| male | 1.05 | [1.03-1.07] | <0.001 |
| married | 0.97 | [0.95-0.99] | <0.001 |
| race_nhb | 0.91 | [0.88-0.94] | <0.001 |
| race_nhapi | 0.79 | [0.77-0.81] | <0.001 |
| race_hispanic | 0.92 | [0.90-0.94] | <0.001 |
| stage_2 | 1.22 | [1.19-1.25] | <0.001 |
| stage_3 | 1.67 | [1.62-1.71] | <0.001 |
| stage_4 | 1.74 | [1.70-1.78] | <0.001 |
| grade_poor | 1.04 | [1.02-1.06] | <0.001 |
| is_icc | 1.29 | [1.26-1.32] | <0.001 |
| surg_local_destruction | 0.26 | [0.25-0.27] | <0.001 |
| surg_segmental_resection | 0.23 | [0.22-0.25] | <0.001 |
| surg_larger_resection | 0.24 | [0.22-0.27] | <0.001 |
| surg_transplant | 0.15 | [0.14-0.17] | <0.001 |
| chemotherapy | 0.57 | [0.56-0.59] | <0.001 |
| radiation | 0.52 | [0.50-0.53] | <0.001 |
| cirrhosis | 0.82 | [0.80-0.85] | <0.001 |

C-index=0.739

## Stratification Analysis

| Strata | NoSurg N | LD HR[CI] | SR HR[CI] | LR HR[CI] | LT HR[CI] |
|--------|---------|-----------|-----------|-----------|-----------|
**age_band** | | | | | |
| 65-69 | 20747 | 0.29[0.27-0.30] | 0.33[0.30-0.36] | 0.35[0.29-0.41] | 0.14[0.13-0.16] |
| 70-74 | 14310 | 0.28[0.27-0.30] | 0.27[0.24-0.31] | 0.25[0.21-0.31] | 0.19[0.16-0.24] |
| 75-79 | 11449 | 0.28[0.26-0.30] | 0.29[0.25-0.33] | 0.33[0.26-0.42] | 0.31[0.24-0.40] |
| 80+ | 13315 | 0.26[0.24-0.29] | 0.26[0.22-0.32] | 0.31[0.22-0.44] | 0.23[0.16-0.34] |
**stage_label** | | | | | |
| Stage I | 25002 | 0.27[0.26-0.29] | 0.26[0.24-0.29] | 0.24[0.20-0.29] | 0.13[0.11-0.15] |
| Stage II | 9625 | 0.36[0.34-0.39] | 0.36[0.31-0.40] | 0.37[0.29-0.48] | 0.13[0.11-0.16] |
| Stage III | 9644 | 0.39[0.35-0.42] | 0.35[0.31-0.40] | 0.37[0.31-0.44] | 0.36[0.30-0.43] |
| Stage IV | 15550 | 0.34[0.31-0.36] | 0.31[0.27-0.36] | 0.28[0.21-0.38] | 0.18[0.15-0.23] |
**cancer_type** | | | | | |
| HCC | 43921 | 0.30[0.29-0.30] | 0.31[0.29-0.33] | 0.31[0.27-0.35] | 0.12[0.11-0.13] |
| ICC | 9553 | 0.23[0.21-0.25] | 0.23[0.20-0.26] | 0.25[0.21-0.30] | 0.25[0.21-0.29] |
| Other | 6347 | 0.26[0.22-0.30] | 0.25[0.20-0.33] | 0.25[0.16-0.39] | 0.21[0.15-0.30] |

## Key Insights

### Median Survival by Surgery Type
  Larger Resection: median OS=24m, mortality=65.4%
  Local Destruction: median OS=28m, mortality=62.3%
  Non-surgery: median OS=4m, mortality=88.3%
  Other: median OS=5m, mortality=85.2%
  Segmental Resection: median OS=30m, mortality=65.4%
  Transplant: median OS=42m, mortality=50.4%

### Age-Surgery Benefit
  65-69: mortality gap = 28.9% (NonSurg 85% vs Surg 56%)
  70-74: mortality gap = 24.9% (NonSurg 88% vs Surg 63%)
  75-79: mortality gap = 19.8% (NonSurg 90% vs Surg 71%)
  80+: mortality gap = 17.4% (NonSurg 92% vs Surg 75%)

### HCC vs ICC
  HCC: median OS=8m, surgery rate=23.5%
  ICC: median OS=5m, surgery rate=18.7%
  Other: median OS=2m, surgery rate=9.1%

### Surgery Benefit by Stage (CSS HR)
  Stage 1: surgery CSS HR=0.26, N=34854
  Stage 2: surgery CSS HR=0.31, N=12816
  Stage 3: surgery CSS HR=0.37, N=11318
  Stage 4: surgery CSS HR=0.34, N=17122

## Figure 2: Forest Plot

✓ Fig2 saved
