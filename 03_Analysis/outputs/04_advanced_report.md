# Advanced Analyses — Hepatobiliary Cancer

## 1. SHAP Feature Analysis

✓ SHAP beeswarm saved

## 2. Competing Risk Analysis

| Cause | N | % of Deaths |
|---|---|---|
| Cancer-Specific Death | 51363 | 81.5% |
| Other-Cause Death | 11625 | 18.5% |
| Alive | 13122 | — |
✓ Competing risk CIF saved

## 3. RMST at 60 Months

| Group | RMST (months) | SE |
|---|---|---|
**Surgery Type** | | |

**Stage** | | |

**Age** | | |

## 4. Age-Surgery Benefit Curve

✓ Age-surgery benefit spline saved

## 5. Novel Insights — Hepatobiliary Cancer

### A. HCC vs ICC Treatment Gap

**HCC**: N=57380, Surgery=23.5%
  None: N=43921, median OS=5m
  Local_Destruction: N=9533, median OS=29m
  Segmental_Resection: N=1815, median OS=32m
  Larger_Resection: N=409, median OS=29m
  Transplant: N=1355, median OS=49m
**ICC**: N=11749, Surgery=18.7%
  None: N=9553, median OS=4m
  Local_Destruction: N=926, median OS=22m
  Segmental_Resection: N=559, median OS=27m
  Larger_Resection: N=251, median OS=20m
  Transplant: N=361, median OS=20m

### B. Cirrhosis Paradox

Cirrhosis=0: surgery rate=20.6%, median OS=6m, mortality=83.8%
Cirrhosis=1: surgery rate=29.5%, median OS=14m, mortality=72.4%

### C. Temporal Trends (2000-2022)

| Year | N | Surgery % | Transplant % | Median OS |
|---|---|---|---|---|
| 2000 | 2089 | 14% | 1.1% | 3m |
| 2001 | 2158 | 14% | 1.2% | 3m |
| 2002 | 2244 | 15% | 2.0% | 3m |
| 2003 | 2184 | 18% | 1.6% | 3m |
| 2004 | 1877 | 23% | 2.8% | 5m |
| 2005 | 2023 | 25% | 3.1% | 5m |
| 2006 | 2057 | 25% | 2.5% | 5m |
| 2007 | 2316 | 22% | 2.0% | 6m |
| 2008 | 2458 | 21% | 2.8% | 6m |
| 2009 | 2591 | 21% | 2.6% | 7m |
| 2010 | 2703 | 21% | 2.1% | 7m |
| 2011 | 2763 | 20% | 2.6% | 8m |
| 2012 | 3099 | 21% | 2.2% | 8m |
| 2013 | 3332 | 23% | 2.1% | 9m |
| 2014 | 3669 | 21% | 2.5% | 10m |
| 2015 | 4026 | 23% | 2.6% | 10m |
| 2016 | 4169 | 23% | 2.0% | 11m |
| 2017 | 4516 | 23% | 2.3% | 11m |
| 2018 | 4738 | 24% | 2.8% | 11m |
| 2019 | 5228 | 24% | 2.4% | 12m |
| 2020 | 4921 | 23% | 2.2% | 10m |
| 2021 | 5588 | 21% | 2.6% | 10m |
| 2022 | 5361 | 20% | 2.5% | 3m |
✓ Temporal trends saved

### D. Segmental vs Larger Resection: HCC vs ICC


**HCC**
  Segmental_Resection: N=1815, median OS=32m, mortality=65.1%, mean age=72.7
  Larger_Resection: N=409, median OS=29m, mortality=63.1%, mean age=72.2
  Segmental_Resection adjusted CSS HR: 0.26 [0.24-0.28]
  Larger_Resection adjusted CSS HR: 0.26 [0.22-0.30]

**ICC**
  Segmental_Resection: N=559, median OS=27m, mortality=64.8%, mean age=72.5
  Larger_Resection: N=251, median OS=20m, mortality=67.7%, mean age=72.5
  Segmental_Resection adjusted CSS HR: 0.22 [0.19-0.24]
  Larger_Resection adjusted CSS HR: 0.26 [0.22-0.31]

### E. Significant Interactions (Cox with interactions)

  surgery_any: HR=0.27, p=<0.001
  age_c: HR=1.02, p=<0.001
  stage_2: HR=1.20, p=<0.001
  stage_3: HR=1.66, p=<0.001
  stage_4: HR=1.87, p=<0.001
  is_icc: HR=1.32, p=<0.001
  grade_poor: HR=0.97, p=0.002
  male: HR=1.03, p=<0.001
  chemotherapy: HR=0.60, p=<0.001
  surg_x_age: HR=1.00, p=0.081
  surg_x_icc: HR=1.00, p=0.916
  surg_x_stage4: HR=1.07, p=0.046

  C-index with interactions: 0.717
