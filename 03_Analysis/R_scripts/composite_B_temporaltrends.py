"""Composite B: Temporal Trends & COVID Impact (2-panel)
Panel A: Annual case counts (bar) + surgery rate (line)  [from Fig15 panels 1-2]
Panel B: Median OS with COVID-19 shading                  [from Fig15 panel 3 + eFigure F]
Deletes: Fig15_COVID.* (old 3-panel version)
"""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from PIL import Image
import os, warnings; warnings.filterwarnings('ignore')

BASE = r'D:\Researching\SEER\hepatobiliary cancer'
os.chdir(BASE)
FIG = '04_Manuscript/figures'
os.makedirs(FIG, exist_ok=True)

df = pd.read_csv(r'02_Data\cleaned\hepatobiliary_elderly_clean.csv')
df['surgery_any'] = (df['surgery_type'] != 'None').astype(int)

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.size': 7,
    'axes.titlesize': 7.5, 'axes.labelsize': 7,
    'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5,
    'legend.fontsize': 6, 'figure.dpi': 300,
    'axes.linewidth': 0.4, 'xtick.major.width': 0.35, 'ytick.major.width': 0.35,
})

def sax(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(0.35)
    ax.spines['bottom'].set_linewidth(0.35)
    ax.tick_params(width=0.35)
    ax.set_facecolor('white')

yrs = np.arange(2004, 2023)
n_vals, surg_vals, os_vals = [], [], []
for y in yrs:
    sub = df[df['year'] == y]
    n_vals.append(len(sub))
    surg_vals.append(sub['surgery_any'].mean() * 100)
    os_vals.append(sub['surv_months'].median())

fig = plt.figure(figsize=(6.85, 3.8))
fig.patch.set_facecolor('white')
gs = GridSpec(1, 2, wspace=0.35, left=0.10, right=0.97, bottom=0.17, top=0.90)

# Panel A: Annual cases (bars, left y) + Surgery rate (line, right y)
ax_a1 = fig.add_subplot(gs[0, 0])
ax_a2 = ax_a1.twinx()
ax_a1.bar(yrs, n_vals, width=0.6, color='#bdc3c7', alpha=0.7, zorder=1, label='Cases')
ax_a1.set_ylabel('Patients (N)', fontsize=7, color='#555')
ax_a1.tick_params(axis='y', labelsize=6.5, colors='#555')
ax_a2.plot(yrs, surg_vals, 'o-', color='#2980b9', lw=1.5, markersize=4, zorder=3, label='Surgery %')
ax_a2.set_ylabel('Surgery Rate (%)', fontsize=7, color='#2980b9')
ax_a2.tick_params(axis='y', labelsize=6.5, colors='#2980b9')
ax_a1.set_xlabel('Year', fontsize=7)
ax_a1.set_title('A. Case Volume & Surgery Rate', fontweight='bold', fontsize=7.5, color='#1a1a1a', loc='left')
lines1, l1 = ax_a1.get_legend_handles_labels(); lines2, l2 = ax_a2.get_legend_handles_labels()
ax_a1.legend(lines1+lines2, l1+l2, frameon=False, fontsize=6)
sax(ax_a1)
ax_a2.spines['top'].set_visible(False)
ax_a2.spines['right'].set_linewidth(0.35)

# Panel B: Median OS with COVID shading
ax_b = fig.add_subplot(gs[0, 1])
ax_b.plot(yrs, os_vals, 's-', color='#c0392b', lw=1.5, markersize=4, zorder=3)
ax_b.axvspan(2019.5, 2022.5, alpha=0.12, color='#c0392b', zorder=0, label='COVID-19')
# Annotations
for yr, val in zip(yrs, os_vals):
    if yr in [2004, 2019, 2020]:
        ax_b.annotate(f'{int(val)}m', (yr, val), textcoords='offset points', xytext=(0, 10),
                      ha='center', fontsize=6, color='#c0392b', fontweight='bold')
ax_b.set_xlabel('Year', fontsize=7); ax_b.set_ylabel('Median OS (months)', fontsize=7)
ax_b.set_title('B. Median OS Trend & COVID-19', fontweight='bold', fontsize=7.5, color='#1a1a1a', loc='left')
ax_b.legend(frameon=False, fontsize=6); sax(ax_b)

name = 'Fig_TemporalTrends_Composite'
fig.savefig(os.path.join(FIG, name+'.png'), dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(os.path.join(FIG, name+'.pdf'), bbox_inches='tight', facecolor='white')
Image.open(os.path.join(FIG, name+'.png')).convert('RGB').save(os.path.join(FIG, name+'.tiff'), 'TIFF', compression='tiff_lzw', dpi=(300,300))
plt.close()

# Delete old Fig15
for old in ['Fig15_COVID.png','Fig15_COVID.pdf','Fig15_COVID.tiff']:
    try: os.remove(os.path.join(FIG, old))
    except: pass
try: os.remove(os.path.join(BASE, '03_Analysis','figures','Fig15_COVID.png'))
except: pass

print(f'  ✓ {name} — replaces Fig15_COVID + eFigure F')
