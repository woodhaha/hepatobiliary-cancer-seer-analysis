"""Composite: AFP (Fig13) + Race (Fig14) + new AFP Pos vs Neg panel
3×3 grid:
  Row 1: AFP Neg (Surg vs No) | AFP Pos (Surg vs No) | AFP Pos vs Neg (overall)
  Row 2: NHW | NHB | NHAPI
  Row 3: Hispanic (center) | blank | blank (hidden)

Each panel: KM curves, median dashed lines, logrank p-value
"""
import pandas as pd, numpy as np
import matplotlib as mpl; mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test
import os, warnings; warnings.filterwarnings('ignore')

BASE = r'D:\Researching\SEER\hepatobiliary cancer'
os.chdir(BASE)
ASO_DIR = '04_Manuscript/figures'
os.makedirs(ASO_DIR, exist_ok=True)

df = pd.read_csv(r'02_Data\cleaned\hepatobiliary_elderly_clean.csv')
df['surgery_any'] = (df['surgery_type'] != 'None').astype(int)

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.size': 8,
    'axes.titlesize': 9, 'axes.labelsize': 8,
    'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5,
    'legend.fontsize': 6.5, 'figure.dpi': 300,
})

KM_COLORS = {'Surgery': '#c0392b', 'No Surgery': '#2980b9',
             'AFP+': '#d35400', 'AFP–': '#2980b9',
             'Positive': '#d35400', 'Negative': '#2980b9'}
SURG_LABEL = {1: 'Surgery', 0: 'No Surgery'}

def sax(ax):
    for s in ['top','right']: ax.spines[s].set_visible(False)
    for s in ['left','bottom']: ax.spines[s].set_linewidth(0.4)
    ax.tick_params(width=0.4); ax.set_facecolor('white')
    ax.set_xlim(0, 60); ax.set_xlabel('Months')
    ax.set_ylabel('Overall Survival')

def add_pval(ax, p):
    txt = 'p < 0.001' if p < 0.001 else f'p = {p:.3f}'
    ax.text(0.97, 0.97, txt, transform=ax.transAxes, fontsize=5.5, fontstyle='italic',
            ha='right', va='top',
            bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='#ccc', lw=0.3, alpha=0.85))

def fit_and_plot(ax, data, group_col, group_vals, group_labels, title):
    """Plot KM curves for groups within data, add median lines and p-value."""
    sax(ax)
    ax.axhline(0.5, color='#ccc', ls=':', lw=0.5, zorder=0)
    pval = None
    for val, lbl in zip(group_vals, group_labels):
        grp = data[data[group_col] == val]
        if len(grp) < 20: continue
        kmf = KaplanMeierFitter()
        kmf.fit(grp['surv_months'], grp['vital_dead'], label=f'{lbl} (n={len(grp)})')
        c = KM_COLORS.get(lbl, '#7f8c8d')
        kmf.plot_survival_function(ax=ax, ci_show=False, lw=1.5, color=c)
        med = kmf.median_survival_time_
        if not np.isnan(med):
            ax.axvline(med, color=c, ls='--', lw=0.5, alpha=0.4)
    # Logrank test
    lr = multivariate_logrank_test(data['surv_months'], data[group_col], data['vital_dead'])
    pval = lr.p_value
    add_pval(ax, pval)
    ax.set_title(title, fontweight='bold', fontsize=7.5, loc='left')
    # Legend bottom left outside curves
    if ax.get_legend_handles_labels()[0]:
        ax.legend(frameon=False, fontsize=5, loc='best')

# ── Figure: 3×3 ──────────────────────────────────────────────
fig = plt.figure(figsize=(7.2, 7.5))
fig.patch.set_facecolor('white')
gs = GridSpec(3, 3, hspace=0.32, wspace=0.28, left=0.08, right=0.97, bottom=0.07, top=0.95)

# ── Row 1: AFP panels ─────────────────────────────────────────
# Panel A: AFP Negative
sub = df[df['afp'] == 'Negative']
ax = fig.add_subplot(gs[0, 0])
fit_and_plot(ax, sub, 'surgery_any', [0, 1], ['No Surgery', 'Surgery'],
             f'AFP Negative (n={len(sub)})')

# Panel B: AFP Positive
sub = df[df['afp'] == 'Positive']
ax = fig.add_subplot(gs[0, 1])
fit_and_plot(ax, sub, 'surgery_any', [0, 1], ['No Surgery', 'Surgery'],
             f'AFP Positive (n={len(sub)})')

# Panel C: AFP Pos vs Neg (overall, no surgery split)
sub = df[df['afp'].isin(['Positive', 'Negative'])]
ax = fig.add_subplot(gs[0, 2])
fit_and_plot(ax, sub, 'afp', ['Negative', 'Positive'], ['AFP–', 'AFP+'],
             f'AFP Positive vs Negative (n={len(sub)})')

# ── Row 2: Race panels (NHW, NHB, NHAPI) ─────────────────────
for i, race in enumerate(['NHW', 'NHB', 'NHAPI']):
    sub = df[df['race'] == race]
    ax = fig.add_subplot(gs[1, i])
    fit_and_plot(ax, sub, 'surgery_any', [0, 1], ['No Surgery', 'Surgery'],
                 f'{race} (n={len(sub)})')

# ── Row 3: Hispanic (center), hide side panels ────────────────
sub = df[df['race'] == 'Hispanic']
ax = fig.add_subplot(gs[2, 1])
fit_and_plot(ax, sub, 'surgery_any', [0, 1], ['No Surgery', 'Surgery'],
             f'Hispanic (n={len(sub)})')

for r, c in [(2, 0), (2, 2)]:
    ax = fig.add_subplot(gs[r, c])
    ax.set_visible(False)

# ── Save ──────────────────────────────────────────────────────
name = 'eFig_AFP_Race_Composite'
plt.rcParams['tiff.compression'] = 'tiff_lzw'
fig.savefig(os.path.join(ASO_DIR, name+'.tiff'), dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(os.path.join(ASO_DIR, name+'.pdf'), bbox_inches='tight', facecolor='white')
plt.close()
print(f'  ✓ {name}')
