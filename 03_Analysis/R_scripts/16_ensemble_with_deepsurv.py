"""Add DeepSurv to the ensemble and compute actual 4-model C-index"""
import pandas as pd, numpy as np, torch, torch.nn as nn
import os, warnings; warnings.filterwarnings('ignore')
os.chdir(r'D:\Researching\SEER\hepatobiliary cancer')

df = pd.read_csv(r'02_Data\cleaned\hepatobiliary_elderly_clean.csv')
df['surgery_type'] = df['surgery_type'].fillna('None')
df['surgery_any'] = (df['surgery_type'] != 'None').astype(int)
df['tumor_size'] = df['tumor_size'].fillna(df['tumor_size'].median())

features = ['age_c','male','married','race_nhb','race_nhapi','race_hispanic',
    'stage_2','stage_3','stage_4','grade_poor','is_icc',
    'chemotherapy','radiation','cirrhosis','income_10k','tumor_size','surgery_any']

from sksurv.util import Surv
from sksurv.metrics import concordance_index_censored
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold
import xgboost as xgb
from sksurv.linear_model import CoxPHSurvivalAnalysis
from sksurv.ensemble import RandomSurvivalForest

ml = df[['surv_months','vital_dead']+features].dropna()
X = ml[features].values.astype(np.float32)
y = Surv.from_arrays(ml['vital_dead'].values.astype(bool), ml['surv_months'].values.astype(np.float64))
scaler = StandardScaler(); X_s = scaler.fit_transform(X)

# DeepSurv model
class DeepSurvV1(nn.Module):
    def __init__(self, n_feats):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(n_feats, 64), nn.BatchNorm1d(64), nn.ReLU(), nn.Dropout(0.3),
            nn.Linear(64, 32), nn.BatchNorm1d(32), nn.ReLU(), nn.Dropout(0.3),
            nn.Linear(32, 16), nn.ReLU(), nn.Linear(16, 1))
    def forward(self, x): return self.net(x)

# 5-fold CV with ALL 4 models
kf = KFold(n_splits=5, shuffle=True, random_state=42)
cv_3model = []; cv_4model = []

for fold, (tr, te) in enumerate(kf.split(X_s)):
    X_tr, X_te = X_s[tr], X_s[te]
    y_tr, y_te = y[tr], y[te]

    # Cox
    cox = CoxPHSurvivalAnalysis(alpha=0.01).fit(X_tr, y_tr)
    p_cox = cox.predict(X_te)

    # RSF
    rsf = RandomSurvivalForest(n_estimators=100, min_samples_split=50, min_samples_leaf=20, n_jobs=-1, random_state=fold).fit(X_tr, y_tr)
    p_rsf = rsf.predict(X_te)

    # XGBoost
    dtrain = xgb.DMatrix(X_tr, label=y_tr['time'])
    dtrain.set_float_info('label_lower_bound', y_tr['time'])
    dtrain.set_float_info('label_upper_bound', np.where(y_tr['event'], y_tr['time'], np.inf))
    xgb_m = xgb.train({'objective':'survival:cox','eval_metric':'cox-nloglik',
        'max_depth':4,'learning_rate':0.05,'tree_method':'hist','seed':fold}, dtrain, num_boost_round=150)
    p_xgb = xgb_m.predict(xgb.DMatrix(X_te))

    # DeepSurv (quick train per fold)
    X_tr_t = torch.FloatTensor(X_tr); X_te_t = torch.FloatTensor(X_te)
    y_t = torch.FloatTensor(y_tr['time'].astype(np.float32))
    y_e = torch.FloatTensor(y_tr['event'].astype(np.float32))

    ds = DeepSurvV1(len(features))
    opt = torch.optim.Adam(ds.parameters(), lr=1e-3, weight_decay=1e-5)
    best_loss = float('inf'); best_state = None; no_imp = 0
    for ep in range(200):
        ds.train(); opt.zero_grad()
        risk = ds(X_tr_t)
        # Cox loss with sort by TIME
        r = risk.squeeze(); _, idx = torch.sort(y_t, descending=True)
        r_s, e_s = r[idx], y_e[idx]
        log_cs = torch.logcumsumexp(r_s, dim=0)
        loss = -(r_s[e_s==1] - log_cs[e_s==1]).sum() / max(e_s.sum(), 1)
        loss.backward(); opt.step()
        if loss.item() < best_loss:
            best_loss = loss.item(); best_state = ds.state_dict().copy(); no_imp = 0
        else: no_imp += 1
        if no_imp >= 15: break

    ds.load_state_dict(best_state); ds.eval()
    with torch.no_grad(): p_ds = ds(X_te_t).cpu().squeeze().numpy()

    # Standardize predictions for ensemble
    def stdize(p): return (p - p.mean()) / (p.std() + 1e-8)
    p3 = np.column_stack([stdize(p_cox), stdize(p_rsf), stdize(p_xgb)])
    p4 = np.column_stack([stdize(p_cox), stdize(p_rsf), stdize(p_xgb), stdize(p_ds)])

    c3 = concordance_index_censored(y_te['event'], y_te['time'], p3.mean(axis=1))[0]
    c4 = concordance_index_censored(y_te['event'], y_te['time'], p4.mean(axis=1))[0]
    cv_3model.append(c3); cv_4model.append(c4)

    # Per-model C-indices
    c_cox = concordance_index_censored(y_te['event'], y_te['time'], p_cox)[0]
    c_rsf = concordance_index_censored(y_te['event'], y_te['time'], p_rsf)[0]
    c_xgb = concordance_index_censored(y_te['event'], y_te['time'], p_xgb)[0]
    c_ds = concordance_index_censored(y_te['event'], y_te['time'], p_ds)[0]
    print(f"Fold {fold+1}: Cox={c_cox:.4f} RSF={c_rsf:.4f} XGB={c_xgb:.4f} DeepSurv={c_ds:.4f} | 3-Ens={c3:.4f} 4-Ens={c4:.4f}")

print(f"\n{'='*60}")
print(f"3-model Ensemble (Cox+RSF+XGB): {np.mean(cv_3model):.4f} ± {np.std(cv_3model):.4f}")
print(f"4-model Ensemble (+DeepSurv):   {np.mean(cv_4model):.4f} ± {np.std(cv_4model):.4f}")
print(f"Improvement: {np.mean(cv_4model)-np.mean(cv_3model):+.4f}")
print(f"{'='*60}")
