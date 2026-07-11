# Key Decisions Log — Hepatobiliary Cancer SEER Analysis

> Generated: 2026-07-11 · Project: Surgical Resection in Elderly Hepatobiliary Cancer

## Study Design Decisions

| Date | Decision | Rationale |
|------|----------|-----------|
| 2026-06-26 | Age cut-off ≥65 for "elderly" | Consistent with Zhang et al. (2020) and WHO definition |
| 2026-06-26 | Unified AJCC staging algorithm | Account for 3rd/6th/7th/EOD 2018 staging changes over study period |
| 2026-06-26 | 5-category surgery classification | None/Local/Segmental/Larger/Transplant following SEER codes |
| 2026-06-26 | Temporal split 2004-2017 train / 2018-2022 test | Clinical temporal validation, avoids data leakage |
| 2026-06-26 | PSM 1:1 nearest-neighbor, caliper=0.05 | Standard in surgical SEER literature |
| 2026-06-26 | E-value on PSM HR (0.35) as primary | More conservative than multivariate HR (0.25) per panel review |

## Analysis Decisions

| Date | Decision | Rationale |
|------|----------|-----------|
| 2026-06-26 | DeepSurv fix: sort by time, not risk | Correct Cox partial likelihood implementation |
| 2026-06-26 | 4-model ensemble (Cox+RSF+XGB+DeepSurv) | After DeepSurv fix, includes best-performing model |
| 2026-06-26 | "Cause-specific hazard" not "Fine-Gray" | Code used CoxPHFitter on filtered subsets, not true subdistribution |
| 2026-06-26 | E-value 5.2 (PSM) reported as primary | Panel review C3: multivariate E-value 7.5 was overstated |
| 2026-06-26 | "Predictor range restriction" not "Surgery Feature Paradox" | Panel review m5: more accurate terminology |
| 2026-07-11 | Manuscript word count corrected to ~2,500 | Audit found overstatement (claimed 3,800) |

## Reporting Decisions

| Date | Decision | Rationale |
|------|----------|-----------|
| 2026-06-26 | TRIPOD Type 3 as primary reporting guideline | Prediction model development + external validation |
| 2026-06-26 | STROBE completed as secondary | Cohort study design also requires STROBE |
| 2026-06-26 | ICGC near-chance performance honestly reported | 95% CI crosses 0.5; population drift acknowledged |

## Manuscript Decisions

| Date | Decision | Rationale |
|------|----------|-----------|
| 2026-07-11 | Ref 8 corrected: J Gastrointest Oncol, not J Gastrointest Surg | Audit found incorrect journal name and year |
| 2026-07-11 | Ref 7: full author list + DOI added | Was incomplete |
| 2026-07-11 | Refs [11,12,17,18,21,22] confirmed cited via multi-ref brackets | False positive from simple regex check |
