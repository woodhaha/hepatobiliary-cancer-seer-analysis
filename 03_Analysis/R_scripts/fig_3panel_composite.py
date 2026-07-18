"""3-panel composite: AFP×Surgery | Surgery×Chemo | Race×Surgery
Each panel: median reference lines + median OS annotations
"""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test as mlt
import os

BASE = r'D:\Researching\SEER\hepatobiliary cancer'
os.chdir(BASE)
ASO_DIR = '04_Manuscript/figures'
os.makedirs(ASO_DIR, exist_ok=True)

df = pd.read_csv(r'02_Data\cleaned\hepatobiliary_elderly_clean.csv', keep_default_na=False)
df['surgery_type'] = df['surgery_type'].replace('', 'None')
df['surgery_any'] = (df['surgery_type'] != 'None').astype(int)

df_race = df[df['race'].isin(['NHW', 'NHB', 'NHAPI', 'Hispanic'])].copy()
df_afp = df[df['afp'].isin(['Positive', 'Negative'])].copy()

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.size': 8,
    'axes.titlesize': 9, 'axes.labelsize': 8,
    'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5,
    'legend.fontsize': 6.5, 'figure.dpi': 300,
})

fig, axes = plt.subplots(1, 3, figsize=(11.5, 4.8))
fig.patch.set_facecolor('white')

XMAX = 65  # extend to show all medians (AFP-/Surg = 63)
KM_COLORS = ['#c0392b', '#e67e22', '#2980b9', '#2ecc71']  # for panel A/B
LS_MAP = {0: '-', 1: '--'}

def fmt_med(m):
    return f'{m:.0f}m' if not np.isnan(m) and m < XMAX else 'NR'

def draw_one_panel(ax, sub_df, groups, title_str):
    """Draw KM curves with median tick x-axis→0.5, median box upper-right."""
    kmf = KaplanMeierFitter()
    meds = {}
    for label, mask, color, ls in groups:
        grp = sub_df[mask]
        if len(grp) == 0: continue
        kmf.fit(grp['surv_months'], grp['vital_dead'], label=f'{label} (n={len(grp)})')
        kmf.plot_survival_function(ax=ax, ci_show=False, lw=1.5, color=color, ls=ls)
        med = kmf.median_survival_time_
        meds[label] = med
        if not np.isnan(med):
            x0 = med if med < XMAX else XMAX - 2
            ax.axvline(x0, ymin=0, ymax=0.5, color=color, ls=ls, lw=0.6, alpha=0.35, zorder=0)
    ax.set_title(title_str, fontweight='bold', fontsize=9)
    ax.set_xlim(0, XMAX)
    ax.legend(frameon=False, loc='lower left', fontsize=6.5)
    med_lines = [f'{k}: {fmt_med(v)}' for k, v in sorted(meds.items())]
    ax.text(0.97, 0.82, 'Median OS\n' + '\n'.join(med_lines), transform=ax.transAxes,
            fontsize=6, ha='right', va='top', fontfamily='monospace',
            bbox=dict(boxstyle='round,pad=0.25', fc='white', ec='#bbb', lw=0.3, alpha=0.9))
    return meds

kmf = KaplanMeierFitter()

# ── Panel A: AFP × Surgery ──
ax = axes[0]
groups_a = [
    ('AFP+ NoSurg', (df_afp['afp']=='Positive') & (df_afp['surgery_any']==0), '#c0392b', '-'),
    ('AFP+ Surg',   (df_afp['afp']=='Positive') & (df_afp['surgery_any']==1), '#e67e22', '--'),
    ('AFP- NoSurg', (df_afp['afp']=='Negative') & (df_afp['surgery_any']==0), '#2980b9', '-'),
    ('AFP- Surg',   (df_afp['afp']=='Negative') & (df_afp['surgery_any']==1), '#2ecc71', '--'),
]
draw_one_panel(ax, df_afp, groups_a, 'AFP Status × Surgery')
ax.set_xlabel('Months'); ax.set_ylabel('Overall Survival')
grp_col = df_afp['afp'] + '_' + df_afp['surgery_any'].map({0:'NS',1:'S'})
pval = mlt(df_afp['surv_months'], grp_col, df_afp['vital_dead']).p_value
p_txt = 'p < 0.001' if pval < 0.001 else f'p = {pval:.3f}'
ax.text(0.97, 0.99, f'Logrank {p_txt}', transform=ax.transAxes,
        fontsize=6.5, fontstyle='italic', ha='right', va='top',
        bbox=dict(boxstyle='round,pad=0.12', fc='white', ec='#ccc', lw=0.3, alpha=0.85))
ax.axhline(0.5, color='#ddd', ls=':', lw=0.5, zorder=0)
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(0.4); ax.spines['bottom'].set_linewidth(0.4)
ax.tick_params(width=0.4)
ax.text(-0.18, 1.02, 'A', transform=ax.transAxes, fontsize=12, fontweight='bold', va='bottom')

# ── Panel B: Surgery × Chemotherapy ──
ax = axes[1]
groups_b = [
    ('NSurg NoChemo', (df['surgery_any']==0) & (df['chemotherapy']==0), '#7f8c8d', '-'),
    ('NSurg Chemo',   (df['surgery_any']==0) & (df['chemotherapy']==1), '#95a5a6', '--'),
    ('Surg NoChemo',  (df['surgery_any']==1) & (df['chemotherapy']==0), '#2980b9', '-'),
    ('Surg Chemo',    (df['surgery_any']==1) & (df['chemotherapy']==1), '#2ecc71', '--'),
]
draw_one_panel(ax, df, groups_b, 'Surgery × Chemotherapy')
ax.set_xlabel('Months'); ax.set_ylabel('')
df['g'] = 'S' + df['surgery_any'].astype(str) + '_C' + df['chemotherapy'].astype(str)
pval = mlt(df['surv_months'], df['g'], df['vital_dead']).p_value
p_txt = 'p < 0.001' if pval < 0.001 else f'p = {pval:.3f}'
ax.text(0.97, 0.99, f'Logrank {p_txt}', transform=ax.transAxes,
        fontsize=6.5, fontstyle='italic', ha='right', va='top',
        bbox=dict(boxstyle='round,pad=0.12', fc='white', ec='#ccc', lw=0.3, alpha=0.85))
ax.axhline(0.5, color='#ddd', ls=':', lw=0.5, zorder=0)
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(0.4); ax.spines['bottom'].set_linewidth(0.4)
ax.tick_params(width=0.4)
ax.text(-0.18, 1.02, 'B', transform=ax.transAxes, fontsize=12, fontweight='bold', va='bottom')

# ── Panel C: Race × Surgery ──
ax = axes[2]
race_colors = {'NHW': '#3498db', 'NHB': '#e74c3c', 'NHAPI': '#2ecc71', 'Hispanic': '#f39c12'}
race_order = ['NHW', 'NHB', 'NHAPI', 'Hispanic']
kmf = KaplanMeierFitter()
meds_c = {}
for race in race_order:
    for surg in [0, 1]:
        mask = (df_race['race']==race) & (df_race['surgery_any']==surg)
        grp = df_race[mask]; clr = race_colors[race]; lss = LS_MAP[surg]
        if len(grp) == 0: continue
        sl = 'Surg' if surg else 'NoS'
        kmf.fit(grp['surv_months'], grp['vital_dead'], label=f'{race} {sl} (n={len(grp)})')
        kmf.plot_survival_function(ax=ax, ci_show=False, lw=1.5, color=clr, ls=lss)
        med = kmf.median_survival_time_
        meds_c[f'{race}-{sl}'] = med
        if not np.isnan(med):
            x0 = med if med < XMAX else XMAX - 2
            ax.axvline(x0, ymin=0, ymax=0.5, color=clr, ls=lss, lw=0.6, alpha=0.35, zorder=0)

df_race['g'] = df_race['race'] + '_' + df_race['surgery_any'].astype(str)
pval_c = mlt(df_race['surv_months'], df_race['g'], df_race['vital_dead']).p_value
ax.set_title('Race × Surgery', fontweight='bold', fontsize=9)
ax.set_xlim(0, XMAX); ax.set_xlabel('Months'); ax.set_ylabel('')
ax.legend(frameon=False, loc='lower left', ncol=1, fontsize=6)
med_lines_c = [f'{k}: {fmt_med(v)}' for k, v in sorted(meds_c.items())]
ax.text(0.97, 0.82, 'Median OS\n' + '\n'.join(med_lines_c), transform=ax.transAxes,
        fontsize=6, ha='right', va='top', fontfamily='monospace',
        bbox=dict(boxstyle='round,pad=0.25', fc='white', ec='#bbb', lw=0.3, alpha=0.9))
p_txt_c = 'p < 0.001' if pval_c < 0.001 else f'p = {pval_c:.3f}'
ax.text(0.97, 0.98, f'Logrank {p_txt_c}', transform=ax.transAxes,
        fontsize=7.5, fontstyle='italic', ha='right', va='top',
        bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='#ccc', lw=0.3, alpha=0.85))
ax.axhline(0.5, color='#ddd', ls=':', lw=0.5, zorder=0)
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(0.4); ax.spines['bottom'].set_linewidth(0.4)
ax.tick_params(width=0.4)
ax.text(-0.18, 1.02, 'C', transform=ax.transAxes, fontsize=12, fontweight='bold', va='bottom')

fig.tight_layout(pad=1.5)

name = 'Fig_3Panel_Composite'
plt.rcParams['tiff.compression'] = 'tiff_lzw'
fig.savefig(os.path.join(ASO_DIR, name+'.tiff'), dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(os.path.join(ASO_DIR, name+'.pdf'), bbox_inches='tight', facecolor='white')
plt.close()
print(f'  ✓ {name} — 3-panel composite with median annotations')
