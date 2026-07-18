"""Merge Fig13 (AFP) + Fig14 (Race) into one composite.
Preserves original comparison layout: each panel = Surgery vs No Surgery KM.
Layout: 2 rows × 3 cols
Row 1: AFP Neg | AFP Pos | NHW
Row 2: NHB | NHAPI | Hispanic
"""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
import os, warnings; warnings.filterwarnings('ignore')

BASE = r'D:\Researching\SEER\hepatobiliary cancer'
os.chdir(BASE)
ASO_DIR = '04_Manuscript/figures'
os.makedirs(ASO_DIR, exist_ok=True)

df = pd.read_csv(r'02_Data\cleaned\hepatobiliary_elderly_clean.csv', keep_default_na=False)
df['surgery_type'] = df['surgery_type'].replace('', 'None')
df['surgery_any'] = (df['surgery_type'] != 'None').astype(int)

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.size': 8,
    'axes.titlesize': 9, 'axes.labelsize': 8,
    'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5,
    'legend.fontsize': 6.5, 'figure.dpi': 300,
})

fig, axes = plt.subplots(2, 3, figsize=(7.2, 5.0))
fig.patch.set_facecolor('white')

panels = [
    # (df_filter, title)
    (df[df['afp'] == 'Negative'], 'AFP Negative'),
    (df[df['afp'] == 'Positive'], 'AFP Positive'),
    (df[df['race'] == 'NHW'],    'NHW'),
    (df[df['race'] == 'NHB'],    'NHB'),
    (df[df['race'] == 'NHAPI'],  'NHAPI'),
    (df[df['race'] == 'Hispanic'], 'Hispanic'),
]

kmf = KaplanMeierFitter()

for ax, (sub, title) in zip(axes.flat, panels):
    n = len(sub)
    for s_lbl, s_grp in sub.groupby('surgery_any'):
        label = f"{'Surgery' if s_lbl==1 else 'No Surgery'}"
        kmf.fit(s_grp['surv_months'], s_grp['vital_dead'], label=label)
        kmf.plot_survival_function(ax=ax, ci_show=False, lw=1.5)
    ax.set_title(f'{title} (n={n})', fontweight='bold')
    ax.set_xlim(0, 60)
    ax.legend(frameon=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(0.4)
    ax.spines['bottom'].set_linewidth(0.4)
    ax.tick_params(width=0.4)

name = 'Fig_AFP_Race_Composite'
plt.rcParams['tiff.compression'] = 'tiff_lzw'
fig.savefig(os.path.join(ASO_DIR, name+'.tiff'), dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(os.path.join(ASO_DIR, name+'.pdf'), bbox_inches='tight', facecolor='white')
plt.close()
print(f'  ✓ {name} — Fig13 (AFP) + Fig14 (Race) merged')
