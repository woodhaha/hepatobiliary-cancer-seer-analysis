"""4-group KM: AFP (Positive/Negative) × Surgery (Yes/No)
One panel with 4 curves, median annotations, logrank p-value
"""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test
import os

BASE = r'D:\Researching\SEER\hepatobiliary cancer'
os.chdir(BASE)
ASO_DIR = '04_Manuscript/figures'
os.makedirs(ASO_DIR, exist_ok=True)

df = pd.read_csv(r'02_Data\cleaned\hepatobiliary_elderly_clean.csv', keep_default_na=False)
df['surgery_type'] = df['surgery_type'].replace('', 'None')
df['surgery_any'] = (df['surgery_type'] != 'None').astype(int)

# Filter to known AFP
df = df[df['afp'].isin(['Positive', 'Negative'])].copy()

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.size': 8,
    'axes.titlesize': 9, 'axes.labelsize': 8,
    'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5,
    'legend.fontsize': 6.5, 'figure.dpi': 300,
})

fig, ax = plt.subplots(figsize=(6.5, 5.0))
fig.patch.set_facecolor('white')

# Define 4 groups
groups = [
    ('AFP+ / No Surgery',     (df['afp']=='Positive') & (df['surgery_any']==0), '#c0392b', '-'),
    ('AFP+ / Surgery',        (df['afp']=='Positive') & (df['surgery_any']==1), '#e67e22', '--'),
    ('AFP− / No Surgery',     (df['afp']=='Negative') & (df['surgery_any']==0), '#2980b9', '-'),
    ('AFP− / Surgery',        (df['afp']=='Negative') & (df['surgery_any']==1), '#2ecc71', '--'),
]

kmf = KaplanMeierFitter()

for label, mask, color, ls in groups:
    grp = df[mask]
    if len(grp) == 0: continue
    kmf.fit(grp['surv_months'], grp['vital_dead'], label=f'{label} (n={len(grp)})')
    kmf.plot_survival_function(ax=ax, ci_show=False, lw=1.5, color=color, ls=ls)
    med = kmf.median_survival_time_
    if not np.isnan(med):
        ax.axvline(med, color=color, ls=ls, lw=0.6, alpha=0.4)

# Logrank across all 4 groups
from lifelines.statistics import multivariate_logrank_test as mlt
df['group'] = df['afp'] + '_' + df['surgery_any'].map({0:'NoSurg',1:'Surg'})
pval = mlt(df['surv_months'], df['group'], df['vital_dead']).p_value

ax.set_xlim(0, 60)
ax.set_xlabel('Months')
ax.set_ylabel('Overall Survival')
ax.set_title('AFP Status × Surgical Treatment', fontweight='bold', loc='left')
ax.legend(frameon=False, loc='lower left')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(0.4)
ax.spines['bottom'].set_linewidth(0.4)
ax.tick_params(width=0.4)

# Add p-value
p_txt = 'p < 0.001' if pval < 0.001 else f'p = {pval:.3f}'
ax.text(0.97, 0.97, f'Logrank {p_txt}', transform=ax.transAxes,
        fontsize=7, fontstyle='italic', ha='right', va='top',
        bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='#ccc', lw=0.3, alpha=0.85))

# Reference 0.5 line
ax.axhline(0.5, color='#ddd', ls=':', lw=0.5, zorder=0)

name = 'Fig_AFP_Surgery_4Group'
plt.rcParams['tiff.compression'] = 'tiff_lzw'
fig.savefig(os.path.join(ASO_DIR, name+'.tiff'), dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(os.path.join(ASO_DIR, name+'.pdf'), bbox_inches='tight', facecolor='white')
plt.close()
print(f'  ✓ {name} — 4-group AFP×Surgery KM')
