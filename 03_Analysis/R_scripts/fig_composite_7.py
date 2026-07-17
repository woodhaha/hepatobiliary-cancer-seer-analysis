"""Fig7_Composite: 6-panel — Model analysis + Competing risk (2×3)
Row A-C: LASSO CV | Coef Paths | DCA
Row D-F: Competing Risk CIF by Surgery | Age | Stage
Replaces: Fig_LASSO_DCA.*, FigS2_CompetingRisk.*
"""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from sklearn.linear_model import LassoCV, lasso_path
from sklearn.preprocessing import StandardScaler
from sksurv.util import Surv
from sksurv.ensemble import RandomSurvivalForest
from lifelines import AalenJohansenFitter
from lifelines.statistics import multivariate_logrank_test
from PIL import Image
import os, warnings; warnings.filterwarnings('ignore')

BASE = r'D:\Researching\SEER\hepatobiliary cancer'
os.chdir(BASE)
FIG = '04_Manuscript/figures'
os.makedirs(FIG, exist_ok=True)

df = pd.read_csv(r'02_Data\cleaned\hepatobiliary_elderly_clean.csv')
df['surgery_type'] = df['surgery_type'].fillna('None')
df['surgery_any'] = (df['surgery_type'] != 'None').astype(int)
df['age_band'] = pd.cut(df['age'], bins=[64,70,75,80,100], labels=['65-69','70-74','75-79','80+'])
df['stage_label'] = df['stage'].map({1:'Stage I',2:'Stage II',3:'Stage III',4:'Stage IV'})
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
    'font.family': 'sans-serif', 'font.size': 7,
    'axes.titlesize': 7.5, 'axes.labelsize': 7,
    'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5,
    'legend.fontsize': 5.5, 'figure.dpi': 300,
    'axes.linewidth': 0.4, 'xtick.major.width': 0.35, 'ytick.major.width': 0.35,
})

def sax(ax):
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(0.35); ax.spines['bottom'].set_linewidth(0.35)
    ax.tick_params(width=0.35); ax.set_facecolor('white')

# ── LASSO data ──
alphas_all = np.logspace(-4, 1, 50)
lasso = LassoCV(cv=5, random_state=42, max_iter=5000, alphas=alphas_all)
lasso.fit(X_s, y_all['time'])
alphas_path, coefs_path, _ = lasso_path(X_s, y_all['time'], alphas=alphas_all)
mse_mean = lasso.mse_path_.mean(axis=1); mse_std = lasso.mse_path_.std(axis=1)
asc = np.argsort(lasso.alphas_); alphas_a = lasso.alphas_[asc]; mse_a = mse_mean[asc]; mse_s = mse_std[asc]
path_asc = np.argsort(alphas_path); alp_p = alphas_path[path_asc]; coef_p = coefs_path[:, path_asc]
min_idx = np.argmin(mse_a); opt_a = alphas_a[min_idx]
se_th = mse_a[min_idx] + mse_s[min_idx]; se_idx = np.where(mse_a[:min_idx] <= se_th)[0]
a1se = alphas_a[se_idx[0]] if len(se_idx) > 0 else opt_a

# ── DCA data ──
rsf = RandomSurvivalForest(n_estimators=100, min_samples_split=50, min_samples_leaf=20, n_jobs=-1, random_state=42).fit(X_s, y_all)
rsf_pred = rsf.predict(X_te)
def net_benefit(pred, y_true, y_time, threshold, t_horizon=60):
    high_risk = pred > np.median(pred); event = (y_time <= t_horizon) & y_true.astype(bool)
    tp = (high_risk & event).sum(); fp = (high_risk & ~event).sum(); n = len(pred)
    return tp/n - fp/n * (threshold/(1-threshold)), event.sum()/n - (~event).sum()/n * (threshold/(1-threshold)), 0
thresholds = np.linspace(0.01, 0.5, 50)
nb_m, nb_a, nb_n = [], [], []
for t in thresholds:
    m, a, n_ = net_benefit(rsf_pred, y_te['event'], y_te['time'], t)
    nb_m.append(m); nb_a.append(a); nb_n.append(0)

# ── Competing risk data ──
df['comp_status'] = 0
df.loc[df['css_dead']==1, 'comp_status'] = 1
df.loc[(df['vital_dead']==1) & (df['css_dead']==0), 'comp_status'] = 2

# ── 2×3 Figure ──
fig = plt.figure(figsize=(6.85, 6.5))
fig.patch.set_facecolor('white')
gs = GridSpec(2, 3, hspace=0.45, wspace=0.35, left=0.08, right=0.97, bottom=0.08, top=0.94)

# ── Row 0: LASSO CV | Coef. Paths | DCA ──
ax = fig.add_subplot(gs[0, 0])
ax.plot(alphas_a, mse_a, color='#2c3e50', lw=1.2, label='CV MSE', zorder=2)
ax.fill_between(alphas_a, mse_a - mse_s, mse_a + mse_s, color='#2980b9', alpha=0.15, label='±1 SD', zorder=1)
ax.axvline(opt_a, color='#c0392b', ls='--', lw=0.8, alpha=0.6, zorder=0)
ax.axvline(a1se, color='#e67e22', ls=':', lw=0.8, alpha=0.6, zorder=0)
ax.scatter(opt_a, mse_a[min_idx], s=25, color='#c0392b', zorder=3, marker='D', edgecolors='white', linewidth=0.3)
ax.set_xscale('log'); ax.set_xlabel('Log(λ)', fontsize=6.5); ax.set_ylabel('MSE', fontsize=6.5)
ax.set_title('A. LASSO CV', fontweight='bold', fontsize=7, loc='left')
ax.legend(frameon=False, fontsize=5.5, loc='lower left'); sax(ax)

ax = fig.add_subplot(gs[0, 1])
feat_labels = {'age_c':'Age','male':'Male','married':'Married','race_nhb':'Race: NHB','race_nhapi':'Race: NHAPI','race_hispanic':'Race: Hispanic',
    'stage_2':'Stage II','stage_3':'Stage III','stage_4':'Stage IV','grade_poor':'Poor Grade','is_icc':'ICC',
    'chemotherapy':'Chemotherapy','radiation':'Radiation','cirrhosis':'Cirrhosis','income_10k':'Income ($10k)','tumor_size':'Tumor Size','surgery_any':'Surgery'}
feat_names = [feat_labels.get(c, c.replace('_',' ').title()) for c in features_all]
path_colors = ['#2980b9','#e67e22','#27ae60','#c0392b','#8e44ad','#f39c12','#1abc9c','#d35400','#3498db','#2ecc71','#9b59b6','#e74c3c','#34495e','#16a085','#f1c40f','#7f8c8d','#2c3e50']
for i in range(len(features_all)):
    ax.plot(alp_p, coef_p[i, :], lw=1, alpha=0.8, color=path_colors[i], label=feat_names[i])
ax.axhline(0, color='#333', lw=0.3, alpha=0.5)
ax.set_xscale('log'); ax.set_xlabel('Log(λ)', fontsize=6.5); ax.set_ylabel('Coefficient', fontsize=6.5)
ax.set_title('B. Coef. Paths', fontweight='bold', fontsize=7, loc='left')
ax.legend(frameon=False, fontsize=4.5, loc='upper left', ncol=1, labelspacing=0.25, handlelength=1, handletextpad=0.3); sax(ax)

ax = fig.add_subplot(gs[0, 2])
ax.plot(thresholds, nb_m, label='RSF Model', lw=1.5, color='#2980b9')
ax.plot(thresholds, nb_a, label='Treat All', lw=1.5, color='#e67e22', ls='--')
ax.axhline(0, color='#7f8c8d', ls=':', lw=0.8, label='Treat None')
ax.set_xlabel('Threshold Probability', fontsize=6.5); ax.set_ylabel('Net Benefit', fontsize=6.5)
ax.set_title('C. Decision Curve', fontweight='bold', fontsize=7, loc='left')
ax.legend(frameon=False, fontsize=5.5); ax.set_xlim(0, 0.5); sax(ax)

# ── Row 1: Competing Risk CIF (with median incidence + log-rank P) ──
cr_configs = [
    ('surgery_any', {1:'Surgery', 0:'Non-surgery'}, 'D. Surgery vs Non-surgery', 'css_dead'),
    ('age_band', {}, 'E. By Age', 'css_dead'),
    ('stage_label', {'Stage I':'Stage I','Stage II':'Stage II','Stage III':'Stage III','Stage IV':'Stage IV'}, 'F. By Stage', 'css_dead'),
]
for col, (gb, lbl_map, title, ev_col) in enumerate(cr_configs):
    ax = fig.add_subplot(gs[1, col])
    dur_cr, ev_cr, grp_cr = [], [], []
    for name, grp in df.groupby(gb):
        if len(grp) < 100: continue
        label = lbl_map.get(name, str(name))
        ajf = AalenJohansenFitter()
        ajf.fit(grp['surv_months'], grp['comp_status'], event_of_interest=1, label=label)
        cd = ajf.cumulative_density_
        ci_vals = cd.iloc[:, 0].values
        times = cd.index.values
        # Find time to reach CIF = 0.5
        idx_50 = np.searchsorted(ci_vals, 0.5)
        k50 = times[idx_50] if idx_50 < len(times) else None
        label_show = label
        if k50 is not None:
            label_show += f', {int(k50)}mo' if np.isfinite(k50) else ''
        ajf.plot(ax=ax, ci_show=False, lw=1.2, label=label_show)
        if k50 is not None:
            ax.plot([k50, k50], [0, 0.5], '--', lw=0.5, alpha=0.4, zorder=0)
        dur_cr.extend(grp['surv_months'].tolist())
        ev_cr.extend(grp[ev_col].tolist())
        grp_cr.extend([label]*len(grp))
    ax.axhline(0.5, color='#999', ls=':', lw=0.35, alpha=0.5, zorder=0)
    if len(set(grp_cr)) >= 2:
        lr = multivariate_logrank_test(dur_cr, ev_cr, grp_cr)
        pv = lr.p_value
        ax.text(0.97, 0.97, f'Log-rank P {"<0.001" if pv < 0.001 else f"={pv:.3f}"}',
                transform=ax.transAxes, fontsize=6, va='top', ha='right', color='#333', style='italic')
    ax.set_title(title, fontweight='bold', fontsize=7, loc='left')
    ax.set_ylabel('CSS Cumulative Incidence', fontsize=6.5)
    ax.set_xlabel('Months', fontsize=6.5)
    ax.set_xlim(0, 120); ax.set_xticks(np.arange(0, 121, 24))
    ax.legend(frameon=False, fontsize=5.5); sax(ax)

name = 'Fig7_Composite'
fig.savefig(os.path.join(FIG, name+'.png'), dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(os.path.join(FIG, name+'.pdf'), bbox_inches='tight', facecolor='white')
Image.open(os.path.join(FIG, name+'.png')).convert('RGB').save(os.path.join(FIG, name+'.tiff'), 'TIFF', compression='tiff_lzw', dpi=(300,300))
plt.close()

for old in ['Fig_LASSO_DCA.png','Fig_LASSO_DCA.pdf','Fig_LASSO_DCA.tiff',
            'FigS2_CompetingRisk.png','FigS2_CompetingRisk.pdf','FigS2_CompetingRisk.tiff']:
    try: os.remove(os.path.join(FIG, old))
    except: pass

print(f'  ✓ {name} — 6-panel: Model Analysis + Competing Risk')
