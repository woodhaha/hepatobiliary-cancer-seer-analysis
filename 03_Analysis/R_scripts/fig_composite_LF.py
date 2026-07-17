"""Fig_Composite_LF: Landmark + Frailty — 2×3 KM panels
Top row: Landmark conditional survival at 0, 12, 24 months
Bottom row: Frailty strata (Fit, Pre-frail, Frail) surgery vs non-surgery
Replaces: Fig12_Landmark, Fig22_Frailty
"""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from lifelines import KaplanMeierFitter, CoxPHFitter
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

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.size': 7,
    'axes.titlesize': 7.5, 'axes.labelsize': 7,
    'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5,
    'legend.fontsize': 6, 'figure.dpi': 300,
})

def sax(ax):
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(0.35); ax.spines['bottom'].set_linewidth(0.35)
    ax.tick_params(width=0.35); ax.set_facecolor('white')

kmf = KaplanMeierFitter()

fig = plt.figure(figsize=(6.85, 5.5))
fig.patch.set_facecolor('white')
gs = GridSpec(2, 3, hspace=0.42, wspace=0.30, left=0.08, right=0.97, bottom=0.10, top=0.94)

# ── Top row: Landmark conditional survival ──
landmarks = [0, 12, 24]
for i, lm in enumerate(landmarks):
    ax = fig.add_subplot(gs[0, i])
    survivors = df[df['surv_months'] > lm].copy()
    survivors['surv_from_lm'] = survivors['surv_months'] - lm
    dur_l, ev_l, grp_l = [], [], []
    spl_med = {}
    for label, grp in survivors.groupby('surgery_any'):
        grp = grp.dropna(subset=['surv_from_lm','vital_dead'])
        if len(grp) < 10: continue
        kmf.fit(grp['surv_from_lm'], grp['vital_dead'], label='Surgery' if label else 'Non-surgery')
        m = kmf.median_survival_time_
        lbl = f"{'Surgery' if label else 'Non-surgery'} (n={len(grp)})"
        lbl += f', {int(m)}mo' if np.isfinite(m) else ', NR'
        spl_med['Surgery' if label else 'Non-surgery'] = m
        kmf.plot_survival_function(ax=ax, ci_show=True, ci_alpha=0.12, lw=1.2,
                                    color='#2980b9' if label else '#7f8c8d', label=lbl)
        dur_l.extend(grp['surv_from_lm'].tolist())
        ev_l.extend(grp['vital_dead'].tolist())
        grp_l.extend(['Surgery' if label else 'Non-surgery'] * len(grp))
    ax.axhline(0.5, color='#999', ls=':', lw=0.35, alpha=0.5, zorder=0)
    for nm, clr in [('Surgery','#2980b9'),('Non-surgery','#7f8c8d')]:
        m = spl_med.get(nm)
        if m is not None and np.isfinite(m):
            ax.plot([m, m], [0, 0.5], '--', lw=0.5, color=clr, alpha=0.4, zorder=0)
    if len(set(grp_l)) >= 2:
        lr = multivariate_logrank_test(dur_l, ev_l, grp_l)
        pv = lr.p_value
        ax.text(0.98, 0.98, f'Log-rank P {"<0.001" if pv < 0.001 else f"={pv:.3f}"}',
                transform=ax.transAxes, fontsize=6, va='top', ha='right', color='#333', style='italic')
    lbl = f'{chr(65+i)}. Conditional Survival: {lm}mo' if lm > 0 else f'{chr(65+i)}. Conditional Survival: 0mo'
    ax.set_title(lbl, fontweight='bold', fontsize=7.5, loc='left')
    ax.set_xlabel(f'Months from {lm}m', fontsize=6.5)
    ax.set_xlim(0, 60); ax.set_xticks(np.arange(0, 61, 12))
    ax.legend(frameon=False, fontsize=5.5, loc='lower left')
    sax(ax)
df['frailty_proxy'] = (
    (df['age'] >= 80).astype(int) * 2 +
    (df['age'] >= 75).astype(int) * 1 +
    (df['stage'] >= 3).astype(int) * 2 +
    (df['grade_poor']).astype(int) * 1
)
df['frailty_group'] = pd.cut(df['frailty_proxy'], bins=[-1, 1, 3, 10],
                             labels=['Fit','Pre-frail','Frail'])

for i, fg in enumerate(['Fit','Pre-frail','Frail']):
    ax = fig.add_subplot(gs[1, i])
    sub = df[df['frailty_group'] == fg]
    dur_f, ev_f, grp_f = [], [], []
    f_med = {}
    for s_lbl, s_sub in sub.groupby('surgery_any'):
        s_sub = s_sub.dropna(subset=['surv_months','vital_dead'])
        if len(s_sub) < 10: continue
        kmf.fit(s_sub['surv_months'], s_sub['vital_dead'], label='Surgery' if s_lbl else 'Non-surgery')
        m = kmf.median_survival_time_
        lbl = f"{'Surgery' if s_lbl else 'Non-surgery'} (n={len(s_sub)})"
        lbl += f', {int(m)}mo' if np.isfinite(m) else ', NR'
        f_med['Surgery' if s_lbl else 'Non-surgery'] = m
        kmf.plot_survival_function(ax=ax, ci_show=True, ci_alpha=0.12, lw=1.2,
                                    color='#2980b9' if s_lbl else '#7f8c8d', label=lbl)
        dur_f.extend(s_sub['surv_months'].tolist())
        ev_f.extend(s_sub['vital_dead'].tolist())
        grp_f.extend(['Surgery' if s_lbl else 'Non-surgery'] * len(s_sub))
    ax.axhline(0.5, color='#999', ls=':', lw=0.35, alpha=0.5, zorder=0)
    for nm, clr in [('Surgery','#2980b9'),('Non-surgery','#7f8c8d')]:
        m = f_med.get(nm)
        if m is not None and np.isfinite(m):
            ax.plot([m, m], [0, 0.5], '--', lw=0.5, color=clr, alpha=0.4, zorder=0)
    if len(set(grp_f)) >= 2:
        lr = multivariate_logrank_test(dur_f, ev_f, grp_f)
        pv = lr.p_value
        ax.text(0.98, 0.98, f'Log-rank P {"<0.001" if pv < 0.001 else f"={pv:.3f}"}',
                transform=ax.transAxes, fontsize=6, va='top', ha='right', color='#333', style='italic')
    ax.set_title(f'{chr(68+i)}. Geriatric Risk: {fg}', fontweight='bold', fontsize=7.5, loc='left')
    ax.set_xlabel('Months', fontsize=6.5)
    ax.set_xlim(0, 60); ax.set_xticks(np.arange(0, 61, 12))
    ax.legend(frameon=False, fontsize=5.5, loc='lower left')
    sax(ax)

name = 'Fig_LandmarkFrailty_Composite'
fig.savefig(os.path.join(FIG, name+'.png'), dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(os.path.join(FIG, name+'.pdf'), bbox_inches='tight', facecolor='white')
Image.open(os.path.join(FIG, name+'.png')).convert('RGB').save(os.path.join(FIG, name+'.tiff'), 'TIFF', compression='tiff_lzw', dpi=(300,300))
plt.close()

# Delete old files
for old in ['Fig12_Landmark.png','Fig12_Landmark.pdf','Fig12_Landmark.tiff',
            'Fig22_Frailty.png','Fig22_Frailty.pdf','Fig22_Frailty.tiff']:
    try: os.remove(os.path.join(FIG, old))
    except: pass
try: os.remove(os.path.join(BASE, '03_Analysis','figures','Fig12_Landmark.png'))
except: pass

print(f'  ✓ {name} — merges Landmark + Frailty into 2×3')
