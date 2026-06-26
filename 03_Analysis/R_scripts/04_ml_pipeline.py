"""Phase 6-8: PSM + ML Survival Models + Treatment Recommendation — Hepatobiliary Cancer"""
import pandas as pd
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
import os, warnings, pickle
warnings.filterwarnings('ignore')

os.chdir(r'D:\Researching\SEER\hepatobiliary cancer')
os.makedirs('03_Analysis/figures', exist_ok=True)
os.makedirs('03_Analysis/outputs', exist_ok=True)

df = pd.read_csv(r'02_Data\cleaned\hepatobiliary_elderly_clean.csv')
df['surgery_type'] = df['surgery_type'].fillna('None')
df['surgery_any'] = (df['surgery_type'] != 'None').astype(int)

# ============================================================
# FEATURES
# ============================================================
features = ['age_c','male','married','race_nhb','race_nhapi','race_hispanic',
            'stage_2','stage_3','stage_4','grade_poor','is_icc',
            'chemotherapy','radiation','cirrhosis','income_10k']

# Impute tumor_size
df['tumor_size'] = df['tumor_size'].fillna(df['tumor_size'].median())
features.append('tumor_size')

# Complete cases for ML
ml = df[['surv_months','vital_dead','css_dead','surgery_any','surgery_type','age','stage','cancer_type'] + features].dropna().copy()

print(f"ML dataset: {len(ml)} complete cases")
print(f"Surgery rate: {ml['surgery_any'].mean()*100:.1f}%")

# Train/Test split (temporal: 2004-2017 vs 2018+)
if 'year' in ml.columns:
    train = ml[ml['year'] <= 2017].copy()
    test = ml[ml['year'] >= 2018].copy()
else:
    train, test = train_test_split(ml, test_size=0.3, random_state=42)
print(f"Train: {len(train)}, Test: {len(test)}")

X_train = train[features].values.astype(np.float32)
X_test = test[features].values.astype(np.float32)

# Scale
scaler = StandardScaler()
X_train_s = scaler.fit_transform(X_train)
X_test_s = scaler.transform(X_test)

# Save for later use
np.savez('03_Analysis/outputs/ml_data.npz',
         X_train=X_train, X_test=X_test,
         X_train_s=X_train_s, X_test_s=X_test_s,
         train_surv=train['surv_months'].values,
         train_os_event=train['vital_dead'].values,
         train_css_event=train['css_dead'].values,
         test_surv=test['surv_months'].values,
         test_os_event=test['vital_dead'].values,
         test_css_event=test['css_dead'].values,
         features=features)

with open('03_Analysis/outputs/03_ml_report.md', 'w', encoding='utf-8') as rpt:
    def p(*a, **kw):
        print(*a, **kw)
        print(*a, **kw, file=rpt)

    p("# ML Survival Models — Hepatobiliary Cancer\n")
    p(f"Train: {len(train)}, Test: {len(test)}, Features: {len(features)}\n")

    # ============================================================
    # MODEL 1: Cox PH (baseline)
    # ============================================================
    from sksurv.linear_model import CoxPHSurvivalAnalysis
    from sksurv.ensemble import RandomSurvivalForest, GradientBoostingSurvivalAnalysis
    from sksurv.metrics import concordance_index_censored
    from sksurv.util import Surv

    y_train = Surv.from_arrays(
        event=train['vital_dead'].values.astype(bool),
        time=train['surv_months'].values.astype(np.float64))
    y_test = Surv.from_arrays(
        event=test['vital_dead'].values.astype(bool),
        time=test['surv_months'].values.astype(np.float64))

    results = {}

    # --- COX ---
    p("## Model 1: Cox PH\n")
    cox = CoxPHSurvivalAnalysis(alpha=0.01)
    cox.fit(X_train_s, y_train)
    cox_pred = cox.predict(X_test_s)
    cox_c = concordance_index_censored(y_test['event'], y_test['time'], cox_pred)[0]
    p(f"C-index: {cox_c:.3f}")
    results['Cox PH'] = cox_c

    # --- RSF ---
    p("## Model 2: Random Survival Forest\n")
    rsf = RandomSurvivalForest(n_estimators=100, min_samples_split=50,
                               min_samples_leaf=20, max_features='sqrt',
                               n_jobs=-1, random_state=42)
    rsf.fit(X_train_s, y_train)
    rsf_pred = rsf.predict(X_test_s)
    rsf_c = concordance_index_censored(y_test['event'], y_test['time'], rsf_pred)[0]
    p(f"C-index: {rsf_c:.3f}")
    results['RSF'] = rsf_c

    # --- XGBoost ---
    try:
        import xgboost as xgb
        p("## Model 3: XGBoost (survival:cox)\n")
        dtrain = xgb.DMatrix(X_train_s, label=y_train['time'])
        dtrain.set_float_info('label_lower_bound', y_train['time'])
        dtrain.set_float_info('label_upper_bound',
            np.where(y_train['event'], y_train['time'], np.inf))

        params = {'objective':'survival:cox','eval_metric':'cox-nloglik',
                  'max_depth':4,'learning_rate':0.05,'subsample':0.8,
                  'colsample_bytree':0.8,'min_child_weight':5,'lambda':1.0,
                  'tree_method':'hist','seed':42}
        xgb_model = xgb.train(params, dtrain, num_boost_round=150)

        dtest = xgb.DMatrix(X_test_s)
        xgb_pred = xgb_model.predict(dtest)
        xgb_c = concordance_index_censored(y_test['event'], y_test['time'], xgb_pred)[0]
        p(f"C-index: {xgb_c:.3f}")
        results['XGBoost'] = xgb_c
    except ImportError:
        p("XGBoost not available")
        xgb_model = None
        xgb_c = 0

    # --- Gradient Boosting (very slow, skip for now) ---
    p("## Model 4: Gradient Boosting Survival (skipped — too slow)\n")
    gb_c = 0
    results['GB Survival'] = gb_c

    # --- DeepSurv ---
    try:
        import torch
        import torch.nn as nn

        class DeepSurv(nn.Module):
            def __init__(self, n_feats):
                super().__init__()
                self.net = nn.Sequential(
                    nn.Linear(n_feats, 64), nn.BatchNorm1d(64), nn.ReLU(), nn.Dropout(0.3),
                    nn.Linear(64, 32), nn.BatchNorm1d(32), nn.ReLU(), nn.Dropout(0.3),
                    nn.Linear(32, 16), nn.ReLU(),
                    nn.Linear(16, 1)
                )
            def forward(self, x): return self.net(x)

        def cox_loss(risk, event):
            risk = risk.squeeze()
            _, idx = torch.sort(risk, descending=True)
            risk_sorted = risk[idx]
            event_sorted = event[idx]
            log_cumsum = torch.logcumsumexp(risk_sorted, dim=0)
            loss = -(risk_sorted[event_sorted==1] - log_cumsum[event_sorted==1]).sum()
            return loss / max(event.sum(), 1)

        p("## Model 5: DeepSurv\n")
        device = 'cuda' if torch.cuda.is_available() else 'cpu'
        model = DeepSurv(len(features)).to(device)
        opt = torch.optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-5)

        X_tr = torch.FloatTensor(X_train_s).to(device)
        y_time = torch.FloatTensor(y_train['time'].astype(np.float32)).to(device)
        y_evt = torch.FloatTensor(y_train['event'].astype(np.float32)).to(device)
        X_te = torch.FloatTensor(X_test_s).to(device)

        best_loss = float('inf')
        patience = 15
        no_improve = 0

        for epoch in range(200):
            model.train()
            opt.zero_grad()
            risk = model(X_tr)
            loss = cox_loss(risk, y_evt)
            loss.backward()
            opt.step()

            if loss.item() < best_loss:
                best_loss = loss.item()
                torch.save(model.state_dict(), '03_Analysis/outputs/deepsurv_best.pt')
                no_improve = 0
            else:
                no_improve += 1
            if no_improve >= patience: break

        model.load_state_dict(torch.load('03_Analysis/outputs/deepsurv_best.pt'))
        model.eval()
        with torch.no_grad():
            ds_pred = model(X_te).cpu().squeeze().numpy()
        ds_c = concordance_index_censored(y_test['event'], y_test['time'], ds_pred)[0]
        p(f"C-index: {ds_c:.3f} (epochs: {epoch+1})")
        results['DeepSurv'] = ds_c
    except ImportError:
        p("PyTorch not available")
        results['DeepSurv'] = 0

    # ============================================================
    # MODEL COMPARISON
    # ============================================================
    p("\n## Model Comparison\n")
    p("| Model | C-index |")
    p("|---|---|")
    for name, score in sorted(results.items(), key=lambda x: x[1], reverse=True):
        p(f"| {name} | {score:.3f} |")

    # ============================================================
    # PERMUTATION FEATURE IMPORTANCE (XGBoost or RSF)
    # ============================================================
    p("\n## Feature Importance\n")
    from sklearn.inspection import permutation_importance as perm_imp

    def cindex_score(est, X, y):
        pred = est.predict(X)
        return concordance_index_censored(y['event'], y['time'], pred)[0]

    best_model = rsf  # Default to RSF
    if xgb_c > rsf_c and xgb_model is not None:
        p("Using XGBoost for feature importance")
        # For XGBoost, use a wrapper
        class XGBWrapper:
            def __init__(self, model, n_feats):
                self.model = model
                self.n_feats = n_feats
            def predict(self, X):
                d = xgb.DMatrix(X)
                return self.model.predict(d)
        best_model = XGBWrapper(xgb_model, len(features))
    else:
        p("Using RSF for feature importance")

    try:
        r = perm_imp(best_model, X_test_s[:2000], y_test, scoring=cindex_score,
                     n_repeats=5, random_state=42, n_jobs=-1)
        imp_df = pd.DataFrame({'feature':features,'importance':r.importances_mean,'std':r.importances_std})
        imp_df = imp_df.sort_values('importance', ascending=False)

        p("| Rank | Feature | Importance |")
        p("|---|---|---|")
        for i, row in imp_df.iterrows():
            p(f"| {list(imp_df.index).index(i)+1} | {row['feature']} | {row['importance']:.4f} ± {row['std']:.4f} |")

        # Plot
        fig, ax = plt.subplots(figsize=(10, 6))
        top = imp_df.head(15)
        ax.barh(range(len(top)), top['importance'].values, xerr=top['std'].values, color='steelblue')
        ax.set_yticks(range(len(top)))
        ax.set_yticklabels(top['feature'].values)
        ax.invert_yaxis()
        ax.set_xlabel('Permutation Importance (C-index drop)')
        ax.set_title('Feature Importance — Hepatobiliary Cancer OS', fontweight='bold')
        fig.savefig('03_Analysis/figures/Fig3_FeatureImportance.png', dpi=300, bbox_inches='tight')
        plt.close()
        p("\n✓ Fig3 saved")
    except Exception as e:
        p(f"Feature importance failed: {e}")

    # ============================================================
    # TREATMENT RECOMMENDATION (Counterfactual)
    # ============================================================
    p("\n## Treatment Recommendation\n")
    p("### Counterfactual Analysis (XGBoost/Risk-based)\n")

    # Use RSF to predict risk for each patient under 4 scenarios
    # No Treatment: surgery_any=0, chemotherapy=0, radiation=0
    # Surgery Only: surgery_any=1, surgery_type=dummy, chemo=0, rad=0
    # Chemo Only: surgery_any=0, chemo=1, rad=0
    # Surgery+Chemo: surgery_any=1, chemo=1, rad=0

    # Simplified: modify surgery_any, chemotherapy
    sample_idx = np.random.choice(len(test), min(5000, len(test)), replace=False)
    X_sample = X_test_s[sample_idx].copy()

    scenarios = {
        'No Treatment': {'surgery_any': 0, 'chemotherapy': 0, 'radiation': 0},
        'Surgery Only': {'surgery_any': 1, 'chemotherapy': 0, 'radiation': 0},
        'Chemo Only': {'surgery_any': 0, 'chemotherapy': 1, 'radiation': 0},
        'Surgery+Chemo': {'surgery_any': 1, 'chemotherapy': 1, 'radiation': 0},
    }

    # Get surgery dummies index
    surg_cols = ['surg_local_destruction','surg_segmental_resection',
                 'surg_larger_resection','surg_transplant']

    scenario_preds = {}
    for scenario, changes in scenarios.items():
        X_mod = X_sample.copy()

        # surgery_any
        if 'surgery_any' in features:
            idx = features.index('surgery_any')
            X_mod[:, idx] = changes['surgery_any']

        # chemotherapy
        if 'chemotherapy' in features:
            idx = features.index('chemotherapy')
            X_mod[:, idx] = changes['chemotherapy']

        # radiation
        if 'radiation' in features:
            idx = features.index('radiation')
            X_mod[:, idx] = changes['radiation']

        # For surgery scenarios, set segmental resection as default
        if changes['surgery_any'] == 1:
            for sc in surg_cols:
                if sc in features:
                    idx = features.index(sc)
                    X_mod[:, idx] = 0
            # Default: segmental resection
            if 'surg_segmental_resection' in features:
                X_mod[:, features.index('surg_segmental_resection')] = 1

        scenario_preds[scenario] = rsf.predict(X_mod)

    # Best scenario per patient
    all_risks = np.column_stack(list(scenario_preds.values()))
    best_idx = np.argmin(all_risks, axis=1)
    best_scenario = np.array(list(scenario_preds.keys()))[best_idx]

    p("### Scenario Distribution (N={})\n".format(len(X_sample)))
    p("| Recommended | N | % |")
    p("|---|---|---|")
    for s in scenarios:
        count = (best_scenario == s).sum()
        p(f"| {s} | {count} | {count/len(X_sample)*100:.1f}% |")

    # Surgery benefit by age
    p("\n### Surgery Benefit by Age (Predicted Risk Reduction)\n")
    no_tx = scenario_preds['No Treatment']
    surg = scenario_preds['Surgery Only']
    benefit = no_tx - surg

    ages = test.iloc[sample_idx]['age'].values
    age_bins = pd.cut(ages, bins=[64,70,75,80,100], labels=['65-69','70-74','75-79','80+'])

    p("| Age Band | Mean Benefit | Median Benefit | % Positive |")
    p("|---|---|---|---|")
    for ab in ['65-69','70-74','75-79','80+']:
        mask = age_bins == ab
        b = benefit[mask]
        p(f"| {ab} | {b.mean():.4f} | {np.median(b):.4f} | {(b>0).mean()*100:.1f}% |")

    p("\n## Key ML Insights\n")
    best_model_name = max(results, key=results.get)
    p(f"Best model: {best_model_name} (C-index={results[best_model_name]:.3f})")

    # Risk stratification
    median_risk = np.median(rsf_pred)
    high_risk = rsf_pred > median_risk
    p(f"High risk cutoff: {median_risk:.3f}")
    p(f"High risk survival: {y_test['time'][high_risk].mean():.1f} months")
    p(f"Low risk survival: {y_test['time'][~high_risk].mean():.1f} months")

print("\n=== ML PIPELINE COMPLETE ===")
print("Report: 03_Analysis/outputs/03_ml_report.md")
print("Data: 03_Analysis/outputs/ml_data.npz")
