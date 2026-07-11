# External Validation Report — Hepatobiliary Cancer

## Data Sources
| Cohort | N | HCC% | Age | Male% | Surgery% | Dead% | Med OS |
|---|---|---|---|---|---|---|---|
| SEER Train | 53,277 | 75% | 74 | 65% | 21% | 83% | 7m |
| SEER Test | 22,833 | 75% | 74 | 65% | 22% | 83% | 8m |
| TCGA-LIHC | 269 | 100% | 59 | 70% | 100% | 7% | 22m |
| ICGC-LIRI-JP | 260 | 100% | 67 | 73% | 100% | 32% | 21m |

## Model Performance
| Model | SEER Internal | TCGA-LIHC | ICGC | Δ (SEER-TCGA) |
|---|---|---|---|---|
| XGBoost | 0.686 | 0.592 | 0.547 | +0.094 |
| RSF | 0.685 | 0.567 | 0.551 | +0.117 |
| Cox PH | 0.669 | 0.595 | 0.522 | +0.074 |

## Predictor Range Restriction

