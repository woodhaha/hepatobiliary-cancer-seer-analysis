"""Composite A: Age-Dependent Surgery Benefit (2-panel)
Panel A: RCS-smoothed age × surgery HR (unadjusted + adjusted)  [replaces Fig21]
Panel B: Age × surgery type stratified HR                          [replaces eFigure Panel D]
Deletes: Fig21_RCS_AgeSpline.*, FigS3_AgeSurgeryBenefit.*
"""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.interpolate import CubicSpline
from lifelines import CoxPHFitter
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
    'axes.linewidth': 0.4, 'xtick.major.width': 0.35, 'ytick.major.width': 0.35,
})

def sax(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(0.35)
    ax.spines['bottom'].set_linewidth(0.35)
    ax.tick_params(width=0.35)
    ax.set_facecolor('white')

# ─── Compute age-stratified HRs ───
ages = np.arange(65, 91, 2)
hrs_raw, hrs_adj = [], []
for a in ages:
    sub = df[(df['age'] >= a) & (df['age'] < a + 2)]
    if len(sub) < 200: hrs_raw.append(np.nan); hrs_adj.append(np.nan); continue
    try:
        c = CoxPHFitter(); c.fit(sub[['surv_months','css_dead','surgery_any']], 'surv_months','css_dead')
        hrs_raw.append(np.exp(c.params_['surgery_any']))
        c2 = CoxPHFitter(); c2.fit(sub[['surv_months','css_dead','surgery_any','stage_2','stage_3','stage_4','is_icc','grade_poor','male']],'surv_months','css_dead')
        hrs_adj.append(np.exp(c2.params_['surgery_any']))
    except: hrs_raw.append(np.nan); hrs_adj.append(np.nan)

hrs_raw, hrs_adj = np.array(hrs_raw), np.array(hrs_adj)
valid = ~np.isnan(hrs_adj)
age_fine = np.linspace(65, 89, 200) if valid.sum() > 4 else ages
if valid.sum() > 4:
    cs_raw = CubicSpline(ages[valid], hrs_raw[valid])
    cs_adj = CubicSpline(ages[valid], hrs_adj[valid])

# ─── Age × surgery type HR ───
surg_types = [('Segmental_Resection','Segmental','-'), ('Larger_Resection','Larger','--'), ('Transplant','Transplant','-.')]
surg_hrs, surg_ages = {}, ages
for skey, slab, sls in surg_types:
    h = []
    for a in ages:
        sub = df[(df['age'] >= a) & (df['age'] < a + 2)]
        sub = sub[sub['surgery_type'].isin(['None', skey])]
        if len(sub) < 100: h.append(np.nan); continue
        try:
            c = CoxPHFitter(); c.fit(sub[['surv_months','css_dead','surgery_any']],'surv_months','css_dead')
            h.append(np.exp(c.params_['surgery_any']))
        except: h.append(np.nan)
    surg_hrs[skey] = h

# ─── Figure ───
fig = plt.figure(figsize=(6.85, 3.8))
fig.patch.set_facecolor('white')
gs = GridSpec(1, 2, wspace=0.40, left=0.10, right=0.97, bottom=0.17, top=0.90)

ax = fig.add_subplot(gs[0, 0])
ax.scatter(ages[valid], hrs_raw[valid], s=15, alpha=0.3, color='#7f8c8d', label='Unadjusted', zorder=2)
ax.scatter(ages[valid], hrs_adj[valid], s=20, alpha=0.5, color='#c0392b', label='Adjusted', zorder=2)
if valid.sum() > 4:
    ax.plot(age_fine, cs_raw(age_fine), '--', color='#2980b9', lw=1.5, alpha=0.6)
    ax.plot(age_fine, cs_adj(age_fine), '-', color='#c0392b', lw=2)
ax.axhline(1, color='#333', ls='-', lw=0.3)
ax.set_xlabel('Age', fontsize=7); ax.set_ylabel('Surgery CSS HR', fontsize=7)
ax.set_title('A. Age-Dependent Surgery Benefit', fontweight='bold', fontsize=7.5, color='#1a1a1a', loc='left')
ax.set_ylim(0, 0.8); ax.legend(frameon=False, fontsize=6); sax(ax)
inflection = age_fine[np.argmin(np.gradient(cs_adj(age_fine)))]
ax.text(0.97, 0.97, f'Inflection: age {inflection:.0f}', transform=ax.transAxes, fontsize=6, va='top', ha='right', color='#c0392b', style='italic')

ax = fig.add_subplot(gs[0, 1])
surg_colors = {'Segmental_Resection':'#2980b9','Larger_Resection':'#e67e22','Transplant':'#27ae60'}
for skey, slab, sls in surg_types:
    h = surg_hrs[skey]
    ax.plot(ages, h, 'o-', color=surg_colors[skey], ls=sls, lw=1.5, label=slab, markersize=4)
ax.axhline(1, color='#333', ls='-', lw=0.3)
ax.set_xlabel('Age', fontsize=7); ax.set_ylabel('Surgery CSS HR', fontsize=7)
ax.set_title('B. By Surgery Type', fontweight='bold', fontsize=7.5, color='#1a1a1a', loc='left')
ax.set_ylim(0, 0.8); ax.legend(frameon=False, fontsize=6); sax(ax)

# Save
name = 'Fig_AgeSurgery_Composite'
fig.savefig(os.path.join(FIG, name+'.png'), dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(os.path.join(FIG, name+'.pdf'), bbox_inches='tight', facecolor='white')
Image.open(os.path.join(FIG, name+'.png')).convert('RGB').save(os.path.join(FIG, name+'.tiff'), 'TIFF', compression='tiff_lzw', dpi=(300,300))
plt.close()

# Delete old files
for old in ['Fig21_RCS_AgeSpline.png','Fig21_RCS_AgeSpline.pdf','Fig21_RCS_AgeSpline.tiff','FigS3_AgeSurgeryBenefit.png','FigS3_AgeSurgeryBenefit.pdf','FigS3_AgeSurgeryBenefit.tiff']:
    try: os.remove(os.path.join(FIG, old))
    except: pass
try: os.remove(os.path.join(BASE, '03_Analysis','figures','FigS3_AgeSurgeryBenefit.png'))
except: pass
try: os.remove(os.path.join(BASE, '03_Analysis','figures','FigS3_AgeSurgeryBenefit.pdf'))
except: pass

print(f'  ✓ {name} — replaces Fig21 + FigS3 + eFigure D')
