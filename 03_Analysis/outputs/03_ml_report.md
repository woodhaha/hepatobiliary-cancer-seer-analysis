# ML Survival Models — Hepatobiliary Cancer

Train: 53277, Test: 22833, Features: 16

## Model 1: Cox PH

C-index: 0.669
## Model 2: Random Survival Forest

C-index: 0.685
## Model 3: XGBoost (survival:cox)

C-index: 0.686
## Model 4: Gradient Boosting Survival (skipped — too slow)

## Model 5: DeepSurv

C-index: 0.523 (epochs: 81)

## Model Comparison

| Model | C-index |
|---|---|
| XGBoost | 0.686 |
| RSF | 0.685 |
| Cox PH | 0.669 |
| DeepSurv | 0.523 |
| GB Survival | 0.000 |

## Feature Importance

Using XGBoost for feature importance
Feature importance failed: The 'estimator' parameter of permutation_importance must be an object implementing 'fit'. Got <__main__.XGBWrapper object at 0x0000017A0CAEECF0> instead.

## Treatment Recommendation

### Counterfactual Analysis (XGBoost/Risk-based)

### Scenario Distribution (N=5000)

| Recommended | N | % |
|---|---|---|
| No Treatment | 1515 | 30.3% |
| Surgery Only | 585 | 11.7% |
| Chemo Only | 2148 | 43.0% |
| Surgery+Chemo | 752 | 15.0% |

### Surgery Benefit by Age (Predicted Risk Reduction)

| Age Band | Mean Benefit | Median Benefit | % Positive |
|---|---|---|---|
| 65-69 | -0.0000 | 0.0000 | 25.7% |
| 70-74 | 0.0000 | 0.0000 | 26.5% |
| 75-79 | -0.0000 | 0.0000 | 25.0% |
| 80+ | -0.0000 | 0.0000 | 26.8% |

## Key ML Insights

Best model: XGBoost (C-index=0.686)
High risk cutoff: 532.694
High risk survival: 10.5 months
Low risk survival: 26.6 months
