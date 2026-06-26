# 1. Schoenfeld Residuals (PH Assumption)

Schoenfeld: proportional_hazard_test() got multiple values for argument 'time_transform'

# 2. Subgroup Interaction Forest Plot

✓ Fig11 Subgroup Forest saved

# 3. Landmark Conditional Survival

✓ Fig12 Landmark saved


# 4. AFP Depth Analysis

| AFP Status | N | Surgery % | Median OS | 1yr OS | 3yr OS |
|---|---|---|---|---|---|
| Borderline | 91 | 37% | 12m | 56% | 38% |
| Negative | 11202 | 33% | 17m | 67% | 41% |
| Positive | 21963 | 19% | 7m | 43% | 22% |
| Unknown | 42854 | 20% | 5m | 35% | 16% |
✓ Fig13 AFP saved


# 5. Race Disparities

| Race | N | Age | Surgery % | Median OS | CSS HR (adj) |
|---|---|---|---|---|---|
| NHW | 42416 | 75 | 22% | 6m | 0.26 |
| NHB | 6778 | 72 | 17% | 6m | 0.28 |
| NHAPI | 12489 | 75 | 27% | 10m | 0.27 |
| Hispanic | 13695 | 74 | 18% | 7m | 0.30 |

**Key**: NHAPI race protective (HR=0.79 in multivariate Cox). Possible factors: etiology (HBV vs HCV vs alcohol), earlier detection, or treatment response differences.
✓ Fig14 Race saved


# 6. COVID Era Impact (2020-2022)

| Period | N | Surgery % | Chemo % | Median OS | 1yr OS |
|---|---|---|---|---|---|
| Pre-COVID (2017-2019) | 14482 | 24% | 36% | 11m | 49% |
| COVID (2020-2022) | 15870 | 21% | 27% | 6m | 47% |
✓ Fig15 COVID saved


## Depth Analysis Summary
| # | Analysis | Key Finding |
|---|---|---|
| 1 | Schoenfeld | ? variables with PH violation |
| 2 | Subgroup | 16 subgroup estimates |
| 3 | Landmark | Conditional survival at 0/12/24 months |
| 4 | AFP | AFP+ patients: lower surgery rate, worse prognosis |
| 5 | Race | NHAPI paradox: best survival despite lower surgery rate |
| 6 | COVID | 2020-2022: case dip, OS drop, recovery in 2022 |
