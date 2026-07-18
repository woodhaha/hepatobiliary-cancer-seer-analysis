"""Fig_LASSO_DCA: 3-panel — LASSO CV + Coef Paths + DCA"""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from sklearn.linear_model import LassoCV, lasso_path
from sklearn.preprocessing import StandardScaler
from sksurv.util import Surv
from sksurv.ensemble import RandomSurvivalForest
import os, warnings; warnings.filterwarnings('ignore')

BASE = r'D:\Researching\SEER\hepatobiliary cancer'
os.chdir(BASE)
FIG = '04_Manuscript/figures'
os.makedirs(FIG, exist_ok=True)

df = pd.read_csv(r'02_Data\cleaned\hepatobiliary_elderly_clean.csv')
df['surgery_any'] = (df['surgery_type'] != 'None').astype(int)
features_all = ['age_c','male','married','race_nhb','race_nhapi','race_hispanic',
    'stage_2','stage_3','stage_4','grade_poor','is_icc',
    'chemotherapy','radiation','cirrhosis','income_10k','tumor_size','surgery_any']
ml = df[['surv_months','vital_dead']+features_all].dropna()
y_all = Surv.from_arrays(ml['vital_dead'].values.astype(bool), ml['surv_months'].values.astype(np.float64))
X_all = ml[features_all].values.astype(np.float32)
scaler = StandardScaler(); X_s = scaler.fit_transform(X_all)
n_test = len(ml) // 3
X_tr, X_te = X_s[:-n_test], X_s[-n_test:]
y_tr, y_te = y_all[:-n_test], y_all[-n_test:]

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.size': 8,
    'axes.titlesize': 9, 'axes.labelsize': 8,
    'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5,
    'legend.fontsize': 6.5, 'figure.dpi': 300,
})

def sax(ax):
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(0.4); ax.spines['bottom'].set_linewidth(0.4)
    ax.tick_params(width=0.4); ax.set_facecolor('white')

# ── Panel A: LASSO CV ──
alphas_all = np.logspace(-4, 1, 50)
lasso = LassoCV(cv=5, random_state=42, max_iter=5000, alphas=alphas_all)
lasso.fit(X_s, y_all['time'])
alphas_path, coefs_path, _ = lasso_path(X_s, y_all['time'], alphas=alphas_all)

mse_mean = lasso.mse_path_.mean(axis=1)
mse_std = lasso.mse_path_.std(axis=1)
asc = np.argsort(lasso.alphas_)
alphas_a = lasso.alphas_[asc]
mse_a = mse_mean[asc]
mse_s = mse_std[asc]
path_asc = np.argsort(alphas_path)
alp_p = alphas_path[path_asc]
coef_p = coefs_path[:, path_asc]

min_idx = np.argmin(mse_a)
opt_a = alphas_a[min_idx]
se_th = mse_a[min_idx] + mse_s[min_idx]
se_idx = np.where(mse_a[:min_idx] <= se_th)[0]
a1se = alphas_a[se_idx[0]] if len(se_idx) > 0 else opt_a
opt_pi = np.argmin(np.abs(alp_p - opt_a))
se_pi = np.argmin(np.abs(alp_p - a1se))
nz_min = int((np.abs(coef_p[:, opt_pi]) > 1e-6).sum())
nz_1se = int((np.abs(coef_p[:, se_pi]) > 1e-6).sum())

# ── Panel B: DCA ──
rsf = RandomSurvivalForest(n_estimators=100, min_samples_split=50, min_samples_leaf=20, n_jobs=-1, random_state=42).fit(X_s, y_all)
rsf_pred = rsf.predict(X_te)

def net_benefit(pred, y_true, y_time, threshold, t_horizon=60):
    high_risk = pred > np.median(pred)
    event = (y_time <= t_horizon) & y_true.astype(bool)
    tp = (high_risk & event).sum(); fp = (high_risk & ~event).sum()
    n = len(pred)
    nb_model = tp/n - fp/n * (threshold/(1-threshold))
    nb_all = event.sum()/n - (~event).sum()/n * (threshold/(1-threshold))
    return nb_model, nb_all, 0

thresholds = np.linspace(0.01, 0.5, 50)
nb_m, nb_a, nb_n = [], [], []
for t in thresholds:
    m, a, n_ = net_benefit(rsf_pred, y_te['event'], y_te['time'], t)
    nb_m.append(m); nb_a.append(a); nb_n.append(0)

# ── Figure (1×3: LASSO CV | DCA | Coef Paths) ──
fig = plt.figure(figsize=(6.85, 3.2))
fig.patch.set_facecolor('white')
gs = GridSpec(1, 3, wspace=0.35, left=0.07, right=0.97, bottom=0.20, top=0.90)

# Panel A: LASSO CV
ax = fig.add_subplot(gs[0, 0])
ax.plot(alphas_a, mse_a, color='#2c3e50', lw=1.2, label='CV MSE', zorder=2)
ax.fill_between(alphas_a, mse_a - mse_s, mse_a + mse_s, color='#2980b9', alpha=0.15, label='±1 SD', zorder=1)
ax.axvline(opt_a, color='#c0392b', ls='--', lw=0.8, alpha=0.6, zorder=0)
ax.axvline(a1se, color='#e67e22', ls=':', lw=0.8, alpha=0.6, zorder=0)
ax.scatter(opt_a, mse_a[min_idx], s=25, color='#c0392b', zorder=3, marker='D', edgecolors='white', linewidth=0.3)
ax.set_xscale('log')
ax.set_xlabel('Log(λ)', fontsize=6.5); ax.set_ylabel('MSE', fontsize=6.5)
ax.set_title('A. LASSO CV', fontweight='bold', fontsize=7, loc='left', pad=4)
ax.text(0.97, 0.97, f'λ_min: {opt_a:.4f}', transform=ax.transAxes, fontsize=5, va='top', ha='right', color='#c0392b')
ax.legend(frameon=False, fontsize=5.5, loc='lower left'); sax(ax)

# Panel B: Coefficient Paths
ax = fig.add_subplot(gs[0, 1])
feat_labels = {
    'age_c': 'Age', 'male': 'Male', 'married': 'Married',
    'race_nhb': 'Race: NHB', 'race_nhapi': 'Race: NHAPI', 'race_hispanic': 'Race: Hispanic',
    'stage_2': 'Stage II', 'stage_3': 'Stage III', 'stage_4': 'Stage IV',
    'grade_poor': 'Poor Grade', 'is_icc': 'ICC',
    'chemotherapy': 'Chemotherapy', 'radiation': 'Radiation',
    'cirrhosis': 'Cirrhosis', 'income_10k': 'Income ($10k)', 'tumor_size': 'Tumor Size',
    'surgery_any': 'Surgery',
}
feat_names = [feat_labels.get(c, c.replace('_',' ').title()) for c in features_all]
path_colors = ['#2980b9','#e67e22','#27ae60','#c0392b','#8e44ad','#f39c12','#1abc9c','#d35400','#3498db','#2ecc71','#9b59b6','#e74c3c','#34495e','#16a085','#f1c40f','#7f8c8d','#2c3e50']
for i in range(len(features_all)):
    ax.plot(alp_p, coef_p[i, :], lw=1, alpha=0.8, color=path_colors[i], label=feat_names[i])
ax.axhline(0, color='#333', lw=0.3, alpha=0.5)
ax.set_xscale('log')
ax.set_xlabel('Log(λ)', fontsize=6.5); ax.set_ylabel('Coefficient', fontsize=6.5)
ax.set_title('B. Coef. Paths', fontweight='bold', fontsize=7, loc='left', pad=4)
nz_cnt = [(np.abs(coef_p[:, i]) > 1e-6).sum() for i in range(coef_p.shape[1])]
ax.text(0.97, 0.03, f'Selected: {nz_cnt[-1]}/{len(features_all)}', transform=ax.transAxes, fontsize=5.5, va='bottom', ha='right', color='#555')
ax.legend(frameon=False, fontsize=4.5, loc='upper left', ncol=1, labelspacing=0.25, handlelength=1, handletextpad=0.3)
sax(ax)

# Panel C: DCA
ax = fig.add_subplot(gs[0, 2])
ax.plot(thresholds, nb_m, label='RSF Model', lw=1.5, color='#2980b9')
ax.plot(thresholds, nb_a, label='Treat All', lw=1.5, color='#e67e22', ls='--')
ax.axhline(0, color='#7f8c8d', ls=':', lw=0.8, label='Treat None')
ax.set_xlabel('Threshold Probability', fontsize=6.5); ax.set_ylabel('Net Benefit', fontsize=6.5)
ax.set_title('C. Decision Curve', fontweight='bold', fontsize=7, loc='left', pad=4)
ax.legend(frameon=False, fontsize=5.5); ax.set_xlim(0, 0.5); sax(ax)

name = 'Fig_LASSO_DCA'
plt.rcParams['tiff.compression'] = 'tiff_lzw'
fig.savefig(os.path.join(FIG, name+'.tiff'), dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(os.path.join(FIG, name+'.pdf'), bbox_inches='tight', facecolor='white')
plt.close()

# Delete old files
for old in ['Fig7_LASSO.png','Fig7_LASSO.pdf','Fig7_LASSO.tiff',
            'Fig8_DCA.png','Fig8_DCA.pdf','Fig8_DCA.tiff',
            'FigS_CoefPaths.png','FigS_CoefPaths.pdf','FigS_CoefPaths.tiff']:
    try: os.remove(os.path.join(FIG, old))
    except: pass

print(f'  ✓ {name} — 3-panel: LASSO CV + DCA + Coef Paths')
