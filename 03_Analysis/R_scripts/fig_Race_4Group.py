"""Single-panel KM: 4 races × Surgery/No Surgery = 8 curves
Preserves the Fig14 comparison structure in one panel.
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

df = df[df['race'].isin(['NHW', 'NHB', 'NHAPI', 'Hispanic'])].copy()

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.size': 8,
    'axes.titlesize': 9, 'axes.labelsize': 8,
    'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5,
    'legend.fontsize': 6.5, 'figure.dpi': 300,
})

fig, ax = plt.subplots(figsize=(7.0, 5.5))
fig.patch.set_facecolor('white')

# 8 curves: 4 races × Surgery/No Surgery
race_colors = {
    'NHW': '#3498db',
    'NHB': '#e74c3c',
    'NHAPI': '#2ecc71',
    'Hispanic': '#f39c12',
}
race_order = ['NHW', 'NHB', 'NHAPI', 'Hispanic']

kmf = KaplanMeierFitter()
lines_styles = {0: (0, ()), 1: (0, (4, 2))}  # solid for No Surgery, dashed for Surgery

for race in race_order:
    for surg in [0, 1]:
        mask = (df['race']==race) & (df['surgery_any']==surg)
        grp = df[mask]
        if len(grp) == 0: continue
        surg_label = 'Surgery' if surg else 'No Surgery'
        kmf.fit(grp['surv_months'], grp['vital_dead'], label=f'{race} - {surg_label} (n={len(grp)})')
        kmf.plot_survival_function(ax=ax, ci_show=False, lw=1.5,
                                   color=race_colors[race],
                                   ls=lines_styles[surg])

# Logrank across all groups
df['group'] = df['race'] + '_' + df['surgery_any'].astype(str)
pval = mlt(df['surv_months'], df['group'], df['vital_dead']).p_value

ax.set_xlim(0, 60)
ax.set_xlabel('Months')
ax.set_ylabel('Overall Survival')
ax.set_title('Overall Survival by Race and Surgery Status', fontweight='bold', loc='left')
ax.legend(frameon=False, loc='lower left', ncol=2)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(0.4)
ax.spines['bottom'].set_linewidth(0.4)
ax.tick_params(width=0.4)

p_txt = 'p < 0.001' if pval < 0.001 else f'p = {pval:.3f}'
ax.text(0.97, 0.97, f'Logrank {p_txt}', transform=ax.transAxes,
        fontsize=7, fontstyle='italic', ha='right', va='top',
        bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='#ccc', lw=0.3, alpha=0.85))

ax.axhline(0.5, color='#ddd', ls=':', lw=0.5, zorder=0)

name = 'Fig_Race_4Group'
plt.rcParams['tiff.compression'] = 'tiff_lzw'
fig.savefig(os.path.join(ASO_DIR, name+'.tiff'), dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(os.path.join(ASO_DIR, name+'.pdf'), bbox_inches='tight', facecolor='white')
plt.close()
print(f'  ✓ {name} — 8-curve race×surgery KM')
