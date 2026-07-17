"""Fig_Composite4: Age-Surgery + Temporal Trends — 4-panel
A: RCS age × surgery benefit (adjusted + unadjusted)
B: Age × surgery type stratified HR
C: Annual case volume + surgery rate
D: Median OS with COVID-19 shading
Replaces: Fig_AgeSurgery_Composite, Fig_TemporalTrends_Composite
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
})

def sax(ax):
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(0.35); ax.spines['bottom'].set_linewidth(0.35)
    ax.tick_params(width=0.35); ax.set_facecolor('white')

# ────── Compute Panel A: Age-stratified RCS ──────
ages = np.arange(65, 91, 2)
hrs_raw, hrs_adj = [], []
for a in ages:
    sub = df[(df['age'] >= a) & (df['age'] < a + 2)]
    if len(sub) < 200: hrs_raw.append(np.nan); hrs_adj.append(np.nan); continue
    try:
        c = CoxPHFitter(); c.fit(sub[['surv_months','css_dead','surgery_any']],'surv_months','css_dead')
        hrs_raw.append(np.exp(c.params_['surgery_any']))
        c2 = CoxPHFitter(); c2.fit(sub[['surv_months','css_dead','surgery_any','stage_2','stage_3','stage_4','is_icc','grade_poor','male']],'surv_months','css_dead')
        hrs_adj.append(np.exp(c2.params_['surgery_any']))
    except: hrs_raw.append(np.nan); hrs_adj.append(np.nan)

hrs_raw, hrs_adj = np.array(hrs_raw), np.array(hrs_adj)
valid = ~np.isnan(hrs_adj)
age_fine = np.linspace(65, 89, 200) if valid.sum() > 4 else ages
if valid.sum() > 4:
    cs_raw = CubicSpline(ages[valid], hrs_raw[valid]); cs_adj = CubicSpline(ages[valid], hrs_adj[valid])

# ────── Compute Panel B: Age × surgery type ──────
surg_types = [('Segmental_Resection','Segmental','-'), ('Larger_Resection','Larger','--'), ('Transplant','Transplant','-.')]
surg_hrs = {}
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

# ────── Compute Panel C+D: Temporal data ──────
yrs = np.arange(2004, 2023)
n_vals, surg_vals, os_vals = [], [], []
for y in yrs:
    sub = df[df['year'] == y]
    n_vals.append(len(sub))
    surg_vals.append(sub['surgery_any'].mean() * 100)
    os_vals.append(sub['surv_months'].median())

# ────── Figure ──────
fig = plt.figure(figsize=(6.85, 5.8))
fig.patch.set_facecolor('white')
gs = GridSpec(2, 2, hspace=0.45, wspace=0.40, left=0.10, right=0.97, bottom=0.10, top=0.94)

# Panel A: RCS
ax = fig.add_subplot(gs[0, 0])
ax.scatter(ages[valid], hrs_raw[valid], s=12, alpha=0.25, color='#7f8c8d', label='Unadjusted', zorder=2)
ax.scatter(ages[valid], hrs_adj[valid], s=18, alpha=0.5, color='#c0392b', label='Adjusted', zorder=2)
if valid.sum() > 4:
    ax.plot(age_fine, cs_raw(age_fine), '--', color='#2980b9', lw=1.5, alpha=0.5)
    ax.plot(age_fine, cs_adj(age_fine), '-', color='#c0392b', lw=2)
ax.axhline(1, color='#333', lw=0.3)
ax.set_xlabel('Age', fontsize=7); ax.set_ylabel('Surgery CSS HR', fontsize=7)
ax.set_title('A. Age-Dependent Surgery Benefit', fontweight='bold', fontsize=7.5, loc='left')
ax.set_ylim(0, 0.8); ax.legend(frameon=False, fontsize=6, loc='upper right'); sax(ax)
if valid.sum() > 4:
    infl = age_fine[np.argmin(np.gradient(cs_adj(age_fine)))]
    ax.text(0.97, 0.03, f'Inflection: age {infl:.0f}', transform=ax.transAxes, fontsize=6, va='bottom', ha='right', color='#c0392b', style='italic')

# Panel B: Surgery type
ax = fig.add_subplot(gs[0, 1])
scolor = {'Segmental_Resection':'#2980b9','Larger_Resection':'#e67e22','Transplant':'#27ae60'}
for skey, slab, sls in surg_types:
    ax.plot(ages, surg_hrs[skey], 'o-', color=scolor[skey], ls=sls, lw=1.5, label=slab, markersize=4)
ax.axhline(1, color='#333', lw=0.3)
ax.set_xlabel('Age', fontsize=7); ax.set_ylabel('Surgery CSS HR', fontsize=7)
ax.set_title('B. By Surgery Type', fontweight='bold', fontsize=7.5, loc='left')
ax.set_ylim(0, 0.8); ax.legend(frameon=False, fontsize=6); sax(ax)

# Panel C: Cases + Surgery rate
ax1 = fig.add_subplot(gs[1, 0])
ax2 = ax1.twinx()
ax1.bar(yrs, n_vals, width=0.6, color='#bdc3c7', alpha=0.7, zorder=1, label='Cases')
ax1.set_ylabel('Patients (N)', fontsize=7, color='#555')
ax1.tick_params(axis='y', labelsize=6.5, colors='#555')
ax2.plot(yrs, surg_vals, 'o-', color='#2980b9', lw=1.5, markersize=4, zorder=3, label='Surgery %')
ax2.set_ylabel('Surgery Rate (%)', fontsize=7, color='#2980b9')
ax2.tick_params(axis='y', labelsize=6.5, colors='#2980b9')
ax1.set_xlabel('Year', fontsize=7); ax1.set_xticks(range(2004, 2023, 3))
for _tl in ax1.get_xticklabels(): _tl.set_rotation(35); _tl.set_fontsize(6.5)
ax1.set_title('C. Case Volume & Surgery Rate', fontweight='bold', fontsize=7.5, loc='left')
l1, la = ax1.get_legend_handles_labels(); l2, lb = ax2.get_legend_handles_labels()
ax1.legend(l1+l2, la+lb, frameon=False, fontsize=6, loc='upper right')
sax(ax1); ax2.spines['top'].set_visible(False); ax2.spines['right'].set_linewidth(0.35)

# Panel D: Median OS + COVID
ax = fig.add_subplot(gs[1, 1])
ax.plot(yrs, os_vals, 's-', color='#c0392b', lw=1.5, markersize=4, zorder=3)
ax.axvspan(2019.5, 2022.5, alpha=0.12, color='#c0392b', zorder=0, label='COVID-19')
for yr, val in zip(yrs, os_vals):
    if yr in [2004, 2019, 2020]:
        ax.annotate(f'{int(val)}m', (yr, val), textcoords='offset points', xytext=(0, 10),
                    ha='center', fontsize=6, color='#c0392b', fontweight='bold')
ax.set_xlabel('Year', fontsize=7); ax.set_ylabel('Median OS (months)', fontsize=7)
ax.set_title('D. Median OS & COVID-19', fontweight='bold', fontsize=7.5, loc='left')
ax.set_xticks(range(2004, 2023, 3))
for _tl in ax.get_xticklabels(): _tl.set_rotation(35); _tl.set_fontsize(6.5)
ax.legend(frameon=False, fontsize=6); sax(ax)

# Save
for name in ['Fig_Composite4']:
    fig.savefig(os.path.join(FIG, name+'.png'), dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig(os.path.join(FIG, name+'.pdf'), bbox_inches='tight', facecolor='white')
    Image.open(os.path.join(FIG, name+'.png')).convert('RGB').save(os.path.join(FIG, name+'.tiff'), 'TIFF', compression='tiff_lzw', dpi=(300,300))
plt.close()

# Delete old composites
for old in ['Fig_AgeSurgery_Composite.png','Fig_AgeSurgery_Composite.pdf','Fig_AgeSurgery_Composite.tiff',
            'Fig_TemporalTrends_Composite.png','Fig_TemporalTrends_Composite.pdf','Fig_TemporalTrends_Composite.tiff']:
    try: os.remove(os.path.join(FIG, old))
    except: pass

print('  ✓ Fig_Composite4 — merges AgeSurgery + TemporalTrends')
