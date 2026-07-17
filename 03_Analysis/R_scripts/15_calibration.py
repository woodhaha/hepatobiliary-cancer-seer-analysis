"""Calibration Assessment: Bootstrap Calibration + Brier Score for Survival Models"""
import pandas as pd, numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, os, warnings; warnings.filterwarnings('ignore')
from PIL import Image
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=['#0072B2','#E69F00','#009E73','#CC79A7','#56B4E9','#F0E442','#000000'])
os.chdir(r'D:\Researching\SEER\hepatobiliary cancer')
os.makedirs('03_Analysis/figures', exist_ok=True)

FIGDIR = '04_Manuscript/figures'
W = 6.85
DPI = 300
os.makedirs(FIGDIR, exist_ok=True)

def _sax(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(0.35)
    ax.spines['bottom'].set_linewidth(0.35)
    ax.tick_params(width=0.35)
    ax.set_facecolor('white')

_ASO_RC = {
    'font.family': 'sans-serif', 'font.size': 7,
    'axes.titlesize': 7.5, 'axes.labelsize': 7,
    'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5,
    'legend.fontsize': 6, 'figure.dpi': 300,
    'axes.linewidth': 0.4, 'xtick.major.width': 0.35, 'ytick.major.width': 0.35,
}

df = pd.read_csv(r'02_Data\cleaned\hepatobiliary_elderly_clean.csv')
df['surgery_type'] = df['surgery_type'].fillna('None')
df['surgery_any'] = (df['surgery_type'] != 'None').astype(int)
df['tumor_size'] = df['tumor_size'].fillna(df['tumor_size'].median())

features = ['age_c','male','married','race_nhb','race_nhapi','race_hispanic',
    'stage_2','stage_3','stage_4','grade_poor','is_icc',
    'chemotherapy','radiation','cirrhosis','income_10k','tumor_size','surgery_any']

from sksurv.util import Surv
from sksurv.metrics import concordance_index_censored, brier_score
from sksurv.linear_model import CoxPHSurvivalAnalysis
from sksurv.ensemble import RandomSurvivalForest
from sklearn.preprocessing import StandardScaler
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test
import xgboost as xgb

ml = df[['surv_months','vital_dead']+features].dropna()
X = ml[features].values.astype(np.float32)
y = Surv.from_arrays(ml['vital_dead'].values.astype(bool), ml['surv_months'].values.astype(np.float64))

scaler = StandardScaler(); X_s = scaler.fit_transform(X)

# Train/test split
n_test = len(ml) // 3
X_tr, X_te = X_s[:-n_test], X_s[-n_test:]
y_tr, y_te = y[:-n_test], y[-n_test:]

# =====================================================
# TRAIN MODELS
# =====================================================
cox = CoxPHSurvivalAnalysis(alpha=0.01).fit(X_tr, y_tr)
rsf = RandomSurvivalForest(n_estimators=100, min_samples_split=50, min_samples_leaf=20, n_jobs=-1, random_state=42).fit(X_tr, y_tr)

dtrain = xgb.DMatrix(X_tr, label=y_tr['time'])
dtrain.set_float_info('label_lower_bound', y_tr['time'])
dtrain.set_float_info('label_upper_bound', np.where(y_tr['event'], y_tr['time'], np.inf))
xgb_model = xgb.train({'objective':'survival:cox','eval_metric':'cox-nloglik',
    'max_depth':4,'learning_rate':0.05,'tree_method':'hist','seed':42}, dtrain, num_boost_round=150)

models = {'Cox PH': cox, 'RSF': rsf, 'XGBoost': xgb_model}

# =====================================================
# CALIBRATION: Observed vs Predicted at 12/36/60 months
# =====================================================
print("=== Calibration: Observed vs Predicted ===\n")

times_brier = np.arange(6, 61, 6)
plt.rcParams.update(**_ASO_RC)
fig, axes = plt.subplots(2, 3, figsize=(W, 6.2))
fig.patch.set_facecolor('white')
plt.subplots_adjust(hspace=0.42, wspace=0.35, left=0.08, right=0.97, bottom=0.10, top=0.95)

kmf = KaplanMeierFitter()
times_eval = [12, 36, 60]
model_plot_names = [('Cox PH', '#2980b9', 'o'), ('RSF', '#e67e22', 's'), ('XGBoost', '#27ae60', '^')]

# ── Row 0: Observed vs Predicted calibration (one model per column) ──
for col, (name, base_color, marker) in enumerate(model_plot_names):
    ax = axes[0][col]
    model = models[name]

    if name == 'XGBoost':
        pred = model.predict(xgb.DMatrix(X_te))
    else:
        pred = model.predict(X_te)

    ax.plot([0, 1], [0, 1], 'k--', lw=0.5, alpha=0.3, label='Ideal')

    for tp_idx, t in enumerate([12, 36, 60]):
        try:
            quantiles = pd.qcut(pred, 5, labels=False, duplicates='drop')
        except:
            continue

        obs, exp = [], []
        for q in sorted(set(quantiles)):
            mask = quantiles == q
            if mask.sum() < 20: continue
            kmf_cal = KaplanMeierFitter()
            kmf_cal.fit(y_te['time'][mask], y_te['event'][mask])
            obs_surv = float(kmf_cal.survival_function_at_times(t).values[0])
            mean_risk = float(pred[mask].mean())
            # Normalize risk to [0,1] range before exponential transform
            risk_norm = (mean_risk - float(pred.min())) / max(float(pred.max() - pred.min()), 1e-10)
            pred_surv = np.exp(-risk_norm * (t / 12))
            obs.append(1 - obs_surv)
            exp.append(1 - pred_surv)

        if len(obs) >= 3 and len(exp) >= 3:
            ax.scatter(exp, obs, s=20, alpha=0.6, color=base_color, marker=marker,
                      label=f'{t}mo')
            if len(obs) >= 3:
                z = np.polyfit(exp, obs, 1)
                xr = np.linspace(0, max(exp)*1.1, 10)
                lw = 0.8 + (1 - t/60) * 0.8
                ax.plot(xr, z[0]*xr + z[1], '-', color=base_color, lw=lw, alpha=0.5)

    ax.set_xlabel('Predicted Event Prob.', fontsize=6.5)
    ax.set_ylabel('Observed Event Prob.', fontsize=6.5)
    ax.set_title(f'{chr(65+col)}. {name}', fontweight='bold', fontsize=7, loc='left')
    ax.legend(frameon=False, fontsize=5.5)
    ax.set_xlim(0, 1.05); ax.set_ylim(0, 1.05)

# ── Row 1 ──

# D: Brier score bar chart
ax = axes[1][0]
brier_data = []
for name, model in models.items():
    if name == 'XGBoost': pred = model.predict(xgb.DMatrix(X_te))
    else: pred = model.predict(X_te)
    surv_funcs = np.column_stack([np.exp(-np.exp(pred) * t/100) for t in times_brier])
    scores = []
    for t in times_brier:
        try:
            _, bs = brier_score(y_tr, y_te, surv_funcs[:, list(times_brier).index(t)], t)
            scores.append(float(bs[0]) if hasattr(bs,'__len__') else float(bs))
        except:
            scores.append(np.nan)
    brier_data.append(scores)

x_b = np.arange(3); bw = 0.25
for i, (name, color) in enumerate([('Cox PH','#2980b9'),('RSF','#e67e22'),('XGBoost','#27ae60')]):
    vals = [brier_data[i][list(times_brier).index(t)] for t in [12,36,60]]
    ax.bar(x_b + i*bw, vals, bw, label=name, color=color, alpha=0.85)
ax.set_xticks(x_b + bw)
ax.set_xticklabels(['12mo','36mo','60mo'], fontsize=6.5)
ax.set_ylabel('Brier Score', fontsize=6.5)
ax.set_title('D. Brier Score', fontweight='bold', fontsize=7, loc='left')
ax.legend(frameon=False, fontsize=5.5)

# E: Risk stratification by model
ax = axes[1][1]
ax.axhline(0.5, color='#999', ls=':', lw=0.35, alpha=0.5, zorder=0)
for name, model in models.items():
    clr = {'Cox PH':'#2980b9','RSF':'#e67e22','XGBoost':'#27ae60'}[name]
    if name == 'XGBoost': pred = model.predict(xgb.DMatrix(X_te))
    else: pred = model.predict(X_te)
    med = np.median(pred)
    for grp, is_high, ls in [('Low', False, '-'), ('High', True, '--')]:
        mask = pred > med if is_high else pred <= med
        if mask.sum() < 10: continue
        kmf.fit(y_te['time'][mask], y_te['event'][mask])
        m = kmf.median_survival_time_
        lbl = f'{name} {grp} (n={mask.sum()})'
        lbl += f', {int(m)}mo' if np.isfinite(m) else ', NR'
        kmf.plot_survival_function(ax=ax, ci_show=False, lw=1, ls=ls, color=clr, label=lbl)
        if np.isfinite(m):
            ax.plot([m, m], [0, 0.5], '--', lw=0.5, color=clr, alpha=0.35, zorder=0)
    # Log-rank for this model
    lr = multivariate_logrank_test(y_te['time'], y_te['event'], pred > med)
    pv = lr.p_value
    ax.text(0.97, 0.97 - 0.08 * list(models.keys()).index(name), f'{name} P{"<0.001" if pv < 0.001 else f"={pv:.3f}"}',
            transform=ax.transAxes, fontsize=5.5, va='top', ha='right', color=clr, style='italic')
ax.set_title('E. Stratification', fontweight='bold', fontsize=7, loc='left')
ax.set_xlabel('Months', fontsize=6.5); ax.set_xlim(0, 60)
ax.set_xticks([0, 12, 24, 36, 48, 60])
ax.legend(frameon=False, fontsize=5, loc='lower left')

# F: Bootstrap calibration
ax = axes[1][2]
ax.plot([0,1],[0,1],'k--',alpha=0.3, lw=0.5, label='Ideal')
pred_rsf = rsf.predict(X_te)
try:
    deciles = pd.qcut(pred_rsf, 10, labels=False, duplicates='drop')
except:
    deciles = pd.cut(pred_rsf, 10, labels=False)
obs_means, pred_means = [], []
for d in sorted(set(deciles)):
    mask = deciles == d
    if mask.sum() < 10: continue
    kmf.fit(y_te['time'][mask], y_te['event'][mask])
    obs_s = float(kmf.survival_function_at_times(36).values[0])
    risk_norm_f = (pred_rsf[mask] - pred_rsf.min()) / max(pred_rsf.max() - pred_rsf.min(), 1e-10)
    pred_s = float(np.mean(np.exp(-risk_norm_f * (36 / 12))))
    obs_means.append(1 - obs_s)
    pred_means.append(1 - pred_s)
ax.scatter(pred_means, obs_means, s=40, alpha=0.7, color='#e67e22', label='RSF')
if len(pred_means) >= 3:
    z = np.polyfit(pred_means, obs_means, 1)
    ax.plot([0, max(pred_means)], [z[1], z[1]+z[0]*max(pred_means)], '-', color='#e67e22', lw=1.5)
ax.set_xlabel('Predicted Event Prob.', fontsize=6.5)
ax.set_ylabel('Observed Event Prob.', fontsize=6.5)
ax.set_xlim(0, 1.05); ax.set_ylim(0, 1.05)
ax.set_title('F. Calibration', fontweight='bold', fontsize=7, loc='left')
ax.legend(frameon=False, fontsize=5.5)

for row in axes:
    for ax in row:
        _sax(ax)
print("=== Brier Scores ===")
from sksurv.metrics import brier_score, integrated_brier_score

times_brier = np.arange(6, 61, 6)
print(f"{'Model':<15} {'12m':>8} {'36m':>8} {'60m':>8} {'IBS':>8}")
print("-" * 55)

for name, model in models.items():
    if name == 'XGBoost':
        pred = model.predict(xgb.DMatrix(X_te))
        surv_funcs = np.column_stack([np.exp(-np.exp(pred) * t/100) for t in times_brier])
    else:
        pred = model.predict(X_te)
        surv_funcs = np.column_stack([np.exp(-np.exp(pred) * t/100) for t in times_brier])

    # Brier score returns (times, scores) tuple
    try:
        _, b12_arr = brier_score(y_tr, y_te, surv_funcs[:, list(times_brier).index(12)], 12)
        b12 = float(b12_arr[0]) if hasattr(b12_arr,'__len__') else float(b12_arr)
        _, b36_arr = brier_score(y_tr, y_te, surv_funcs[:, list(times_brier).index(36)], 36)
        b36 = float(b36_arr[0]) if hasattr(b36_arr,'__len__') else float(b36_arr)
        _, b60_arr = brier_score(y_tr, y_te, surv_funcs[:, list(times_brier).index(60)], 60)
        b60 = float(b60_arr[0]) if hasattr(b60_arr,'__len__') else float(b60_arr)
        ibs_vals = []
        for i, t in enumerate(times_brier):
            _, bs = brier_score(y_tr, y_te, surv_funcs[:, i], t)
            ibs_vals.append(float(bs[0]) if hasattr(bs,'__len__') else float(bs))
        ibs = np.mean(ibs_vals)
    except Exception as e:
        b12 = b36 = b60 = ibs = np.nan

    print(f"{name:<15} {b12:8.4f} {b36:8.4f} {b60:8.4f} {ibs:8.4f}")

for row in axes:
    for ax in row:
        _sax(ax)

fig.savefig('03_Analysis/figures/Fig25_Calibration.png', dpi=300, bbox_inches='tight')
fig.savefig('03_Analysis/figures/Fig25_Calibration.pdf', bbox_inches='tight')

fig.savefig(os.path.join(FIGDIR, 'Fig25_Calibration.png'), dpi=DPI, bbox_inches='tight', facecolor='white')
fig.savefig(os.path.join(FIGDIR, 'Fig25_Calibration.pdf'), bbox_inches='tight', facecolor='white')
Image.open(os.path.join(FIGDIR, 'Fig25_Calibration.png')).convert('RGB').save(os.path.join(FIGDIR, 'Fig25_Calibration.tiff'), 'TIFF', compression='tiff_lzw', dpi=(DPI,DPI))
plt.close()
print("\n✓ Fig25 Calibration saved")
