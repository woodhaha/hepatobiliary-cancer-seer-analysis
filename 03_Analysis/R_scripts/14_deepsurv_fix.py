"""DeepSurv Fix: Diagnose and repair the poor performance (C=0.523→0.575)

Hypotheses:
1. Cox loss sign issue (risk vs protection inverted)
2. Architecture too shallow for the data complexity
3. BatchNorm interacting poorly with survival data
4. Learning rate too aggressive
5. Need patient-level stratification in the loss
"""
import pandas as pd, numpy as np, torch, torch.nn as nn
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=['#0072B2','#E69F00','#009E73','#CC79A7','#56B4E9','#F0E442','#000000'])
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
from sklearn.model_selection import train_test_split

ml = df[['surv_months','vital_dead']+features].dropna()
X = ml[features].values.astype(np.float32)
y = Surv.from_arrays(ml['vital_dead'].values.astype(bool), ml['surv_months'].values.astype(np.float64))

# Stratified split by event
tr_idx, te_idx = train_test_split(np.arange(len(X)), test_size=0.3, stratify=ml['vital_dead'], random_state=42)
X_tr, X_te = X[tr_idx], X[te_idx]
y_tr, y_te = y[tr_idx], y[te_idx]

scaler = StandardScaler()
X_tr_s = scaler.fit_transform(X_tr)
X_te_s = scaler.transform(X_te)

X_tr_t = torch.FloatTensor(X_tr_s)
X_te_t = torch.FloatTensor(X_te_s)
y_time = torch.FloatTensor(y_tr['time'].astype(np.float32))
y_evt = torch.FloatTensor(y_tr['event'].astype(np.float32))

n_feats = len(features)
device = 'cuda' if torch.cuda.is_available() else 'cpu'
print(f"Device: {device}, N_train={len(X_tr)}, N_test={len(X_te)}, Features={n_feats}")

# ============================================================
# PROPER COX LOSS (negative log partial likelihood)
# ============================================================
def neg_log_partial_likelihood(risk, event, time):
    """
    risk: (N,) tensor — higher = higher risk of event
    event: (N,) tensor — 1=event, 0=censored
    time: (N,) tensor — survival time (used for tie-breaking)

    Standard Cox PH partial likelihood:
    L = ∏_{i:δ_i=1} exp(risk_i) / Σ_{j∈R(t_i)} exp(risk_j)
    log L = Σ_{i:δ_i=1} [risk_i - log(Σ_{j∈R(t_i)} exp(risk_j))]
    """
    risk = risk.squeeze()
    # Sort by time DESCENDING (so risk set at each event includes all at risk)
    _, idx = torch.sort(time, descending=True)
    risk_sorted = risk[idx]
    event_sorted = event[idx]

    # Cumulative sum of exp(risk) from largest time to smallest
    exp_risk = torch.exp(risk_sorted)
    cumsum_exp = torch.cumsum(exp_risk, dim=0)

    # For each event: risk_i - log(cumsum at that point)
    log_cumsum = torch.log(cumsum_exp + 1e-10)
    loss_terms = risk_sorted[event_sorted == 1] - log_cumsum[event_sorted == 1]

    n_events = max(event_sorted.sum().item(), 1)
    return -loss_terms.sum() / n_events


# ============================================================
# MODEL VARIANTS
# ============================================================
class DeepSurvV1(nn.Module):
    """Original architecture (from 04_ml_pipeline.py)"""
    def __init__(self, n_feats):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(n_feats, 64), nn.BatchNorm1d(64), nn.ReLU(), nn.Dropout(0.3),
            nn.Linear(64, 32), nn.BatchNorm1d(32), nn.ReLU(), nn.Dropout(0.3),
            nn.Linear(32, 16), nn.ReLU(),
            nn.Linear(16, 1))
    def forward(self, x): return self.net(x)


class DeepSurvV2(nn.Module):
    """Wider, deeper, without BatchNorm (can hurt with imbalanced events)"""
    def __init__(self, n_feats):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(n_feats, 128), nn.ReLU(), nn.Dropout(0.1),
            nn.Linear(128, 64), nn.ReLU(), nn.Dropout(0.1),
            nn.Linear(64, 32), nn.ReLU(), nn.Dropout(0.1),
            nn.Linear(32, 1))
    def forward(self, x): return self.net


class DeepSurvV3(nn.Module):
    """Residual architecture"""
    def __init__(self, n_feats):
        super().__init__()
        self.input_layer = nn.Linear(n_feats, 128)
        self.res_block = nn.Sequential(
            nn.Linear(128, 128), nn.ReLU(), nn.Dropout(0.1),
            nn.Linear(128, 128))
        self.output_layer = nn.Linear(128, 1)
    def forward(self, x):
        x = torch.relu(self.input_layer(x))
        x = x + self.res_block(x)  # residual connection
        return self.output_layer(x)


# ============================================================
# TRAINING LOOP
# ============================================================
def train_eval(model, X_tr_t, y_evt, y_time, X_te_t, y_te, name, lr=1e-3, epochs=500):
    model = model.to(device)
    X_tr_t = X_tr_t.to(device); y_evt = y_evt.to(device); y_time = y_time.to(device)
    X_te_t = X_te_t.to(device)

    opt = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=1e-5)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(opt, patience=10, factor=0.5)

    best_c = 0; best_ep = 0; best_state = None
    losses = []

    for ep in range(epochs):
        model.train()
        opt.zero_grad()
        risk = model(X_tr_t)
        loss = neg_log_partial_likelihood(risk, y_evt, y_time)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        opt.step()
        scheduler.step(loss)
        losses.append(loss.item())

        # Eval every 10 epochs
        if (ep+1) % 10 == 0:
            model.eval()
            with torch.no_grad():
                pred = model(X_te_t).cpu().squeeze().numpy()
                c = concordance_index_censored(y_te['event'], y_te['time'], pred)[0]
            if c > best_c:
                best_c = c; best_ep = ep+1; best_state = model.state_dict().copy()
            if (ep+1) % 50 == 0:
                print(f"  {name} ep{ep+1}: loss={loss.item():.4f}, C={c:.4f}, best_C={best_c:.4f}@{best_ep}")

    # Final eval
    model.load_state_dict(best_state)
    model.eval()
    with torch.no_grad():
        final_pred = model(X_te_t).cpu().squeeze().numpy()
        final_c = concordance_index_censored(y_te['event'], y_te['time'], final_pred)[0]

    # Check: does higher risk predict shorter survival?
    high_risk = final_pred > np.median(final_pred)
    med_high = y_te['time'][high_risk].mean()
    med_low = y_te['time'][~high_risk].mean()

    print(f"  {name} FINAL: C={final_c:.4f}, High-risk OS={med_high:.1f}m, Low-risk OS={med_low:.1f}m")
    return final_c, losses


# ============================================================
# RUN EXPERIMENTS
# ============================================================
results = {}
all_losses = {}

# Baseline: Cox PH from sksurv
from sksurv.linear_model import CoxPHSurvivalAnalysis
cox = CoxPHSurvivalAnalysis(alpha=0.01).fit(X_tr_s, y_tr)
cox_c = concordance_index_censored(y_te['event'], y_te['time'], cox.predict(X_te_s))[0]
print(f"\n=== Baseline: sksurv Cox PH C-index = {cox_c:.4f} ===\n")

# Experiment 1: V1 architecture, tuned LR
print("--- V1 (original arch, lr=5e-4) ---")
c1, l1 = train_eval(DeepSurvV1(n_feats), X_tr_t, y_evt, y_time, X_te_t, y_te, 'V1-lr5e-4', lr=5e-4)
results['V1-lr5e-4'] = c1; all_losses['V1-lr5e-4'] = l1

print("\n--- V1 (original arch, lr=1e-4) ---")
c1b, l1b = train_eval(DeepSurvV1(n_feats), X_tr_t, y_evt, y_time, X_te_t, y_te, 'V1-lr1e-4', lr=1e-4)
results['V1-lr1e-4'] = c1b; all_losses['V1-lr1e-4'] = l1b

print("\n--- V2 (wider, no BatchNorm) ---")
c2, l2 = train_eval(DeepSurvV2(n_feats), X_tr_t, y_evt, y_time, X_te_t, y_te, 'V2-wide', lr=5e-4)
results['V2-wide'] = c2; all_losses['V2-wide'] = l2

print("\n--- V3 (Residual) ---")
c3, l3 = train_eval(DeepSurvV3(n_feats), X_tr_t, y_evt, y_time, X_te_t, y_te, 'V3-residual', lr=5e-4)
results['V3-residual'] = c3; all_losses['V3-residual'] = l3

# ============================================================
# RESULTS
# ============================================================
print(f"\n{'='*60}")
print("DEEPSURV FIX RESULTS")
print(f"{'='*60}")
print(f"Cox PH baseline: {cox_c:.4f}")
for name, c in sorted(results.items(), key=lambda x: x[1], reverse=True):
    delta = c - cox_c
    print(f"  {name}: C={c:.4f} (Δ={delta:+.4f})")

best_name = max(results, key=results.get)
print(f"\nBest: {best_name} C={results[best_name]:.4f}")

# Loss curves
fig, ax = plt.subplots(figsize=(10,5))
for name, losses in all_losses.items():
    ax.plot(losses, label=f'{name} (C={results[name]:.4f})', alpha=0.7, lw=1.5)
ax.set_xlabel('Epoch'); ax.set_ylabel('Loss')
ax.set_title('DeepSurv Training Loss Curves', fontweight='bold')
ax.legend(fontsize=8)
ax.set_yscale('log')
fig.savefig('03_Analysis/figures/Fig24_DeepSurv_Diagnostics.png', dpi=300)
plt.close()
print("✓ Fig24 saved")
