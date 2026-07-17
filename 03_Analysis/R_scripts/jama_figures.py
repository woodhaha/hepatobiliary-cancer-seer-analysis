"""JAMA Network Open — Figure Restyling Script
Regenerates Fig1, Fig2, Fig3 with JAMA-compliant styling:
- At-risk tables on KM curves (JAMA requirement)
- Clean white-background minimalist style
- Professional restrained color palette
- TIFF RGB 300dpi output"""

import pandas as pd, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter, CoxPHFitter
import os, warnings; warnings.filterwarnings('ignore')
plt.rcParams.update({
    'font.family': 'sans-serif', 'font.size': 10,
    'axes.titlesize': 11, 'axes.labelsize': 10,
    'xtick.labelsize': 9, 'ytick.labelsize': 9,
    'legend.fontsize': 8, 'figure.dpi': 300,
})

BASE = r'D:\Researching\SEER\hepatobiliary cancer'
os.chdir(BASE)
FIG = os.path.join(BASE, '04_Manuscript', 'figures')
os.makedirs(FIG, exist_ok=True)
os.makedirs('03_Analysis/figures', exist_ok=True)

# JAMA-inspired restrained palette
C = {
    'non_surg': '#7f8c8d',  # gray
    'local':    '#e67e22',  # orange
    'seg':      '#2b7bba',  # JAMA blue
    'larger':   '#8e44ad',  # purple
    'txp':      '#27ae60',  # green
    'hcc':      '#2b7bba',
    'icc':      '#8e44ad',
    'hbi_lo':   '#2b7bba',
    'hbi_mid':  '#e67e22',
    'hbi_hi':   '#c0392b',
    'ref_gray': '#bdc3c7',
}
SURG_COLORS = {
    'None': C['non_surg'], 'Local Destruction': C['local'],
    'Segmental Resection': C['seg'], 'Larger Resection': C['larger'],
    'Transplant': C['txp'],
}

df = pd.read_csv(r'02_Data\cleaned\hepatobiliary_elderly_clean.csv')
df['surgery_type'] = df['surgery_type'].fillna('None').str.replace('_', ' ').replace({'Surgery NOS':'Other','Other Surgery':'Other'})
df['surgery_any'] = (df['surgery_type'] != 'None').astype(int)
df['surgery_vs_none'] = df['surgery_any'].map({1:'Surgery',0:'Non-surgery'})
df['tumor_size'] = df['tumor_size'].fillna(df['tumor_size'].median())
df['age_band'] = pd.cut(df['age'], bins=[64,70,75,80,100], labels=['65-69','70-74','75-79','80+'])
df['stage_label'] = df['stage'].map({1:'I',2:'II',3:'III',4:'IV'}).fillna('?')
# surgery_group: collapse others into 'None' for KM Plot C
df['surgery_group'] = df['surgery_type'].where(df['surgery_type'].isin(
    ['None','Local Destruction','Segmental Resection','Larger Resection','Transplant']), 'Other')
df['survival_months'] = df['surv_months']  # lifelines alias
df['event'] = df['vital_dead']

# ─── helpers ───────────────────────────────────────────────────────
def get_ci(m, col):
    try: return m.confidence_intervals_.loc[col,'95% lower bound'], m.confidence_intervals_.loc[col,'95% upper bound']
    except: return np.exp(m.params_[col]-1.96*m.standard_errors_[col]), np.exp(m.params_[col]+1.96*m.standard_errors_[col])

def save(fig, name):
    path = os.path.join(FIG, name)
    fig.savefig(path + '.png', dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig(path + '.pdf', bbox_inches='tight', facecolor='white')
    # TIFF: ensure RGB mode
    from PIL import Image
    im = Image.open(path + '.png').convert('RGB')
    tiff_path = os.path.join(FIG, name + '.tiff')
    im.save(tiff_path, 'TIFF', compression='tiff_lzw', dpi=(300,300))
    print(f"✓ {name}")

# ─── FIGURE 1: KM with Number-at-Risk Tables ──────────────────────
def add_risk_table(ax, kmfs, times, labels, n0):
    """Add number-at-risk table below a KM plot."""
    rows = []
    for kmf, lb in zip(kmfs, labels):
        max_t = kmf.timeline[-1]
        row = [str(int(n0.get(lb, 0)))]
        for t in times[1:]:
            if t > max_t:
                n = 0
            else:
                # Find at-risk at the closest time >= t
                at_risk_series = kmf.event_table['at_risk']
                n = at_risk_series.loc[at_risk_series.index >= t].iloc[0] if any(at_risk_series.index >= t) else 0
            row.append(str(int(n)))
        rows.append(row)

    ncol = len(times)
    table_ax = ax.figure.add_axes([ax.get_position().x0, ax.get_position().y0-0.14,
                                    ax.get_position().width, 0.10])
    table_ax.axis('off')
    y0 = 0.85; col_w = 1.0 / ncol
    table_ax.text(-0.05, y0, 'No. at risk', fontsize=7, fontweight='bold', va='bottom')
    for j, cl in enumerate([str(int(t)) for t in times]):
        table_ax.text(j*col_w + col_w/2, y0, cl, fontsize=7, ha='center', va='bottom')
    for i, (row, lb) in enumerate(zip(rows, labels)):
        y = y0 - (i+1)*0.20
        table_ax.text(-0.05, y, lb, fontsize=6.5, va='top', color=SURG_COLORS.get(lb, 'black'))
        for j, v in enumerate(row):
            table_ax.text(j*col_w + col_w/2, y, v, fontsize=6.5, ha='center', va='top')

print("Generating Fig1 — Kaplan-Meier with at-risk tables...")
fig, axes = plt.subplots(2, 3, figsize=(14, 12))
fig.patch.set_facecolor('white')

# Collect all kmfits for risk table
all_kmfs = {}
TIMES = [0, 12, 24, 36, 60, 120]

def plot_km(ax, data, group_col, title, n_init, palette):
    kmf = KaplanMeierFitter()
    fits, labs = [], []
    for lb, gr in data.groupby(group_col):
        if len(gr) < 30: continue
        kmf.fit(gr['surv_months'], gr['vital_dead'], label=lb)
        c = palette.get(lb, C['non_surg'])
        kmf.plot_survival_function(ax=ax, ci_show=False, lw=1.5, color=c)
        fits.append(kmf); labs.append(lb)
    ax.set_title(title, fontweight='bold', fontsize=10)
    ax.set_xlabel('Months', fontsize=9); ax.set_ylabel('Survival Probability', fontsize=9)
    ax.set_xlim(0, 120); ax.set_ylim(0, 1.05)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    ax.legend(frameon=False, fontsize=7, loc='lower left')
    return fits, labs

# Panel A
n_init = {lb: len(gr) for lb, gr in df.groupby('surgery_vs_none')}
plot_km(axes[0,0], df, 'surgery_vs_none',
        'A. OS: Surgery vs Non-surgery', n_init, {'Surgery': C['seg'], 'Non-surgery': C['non_surg']})
axes[0,0].legend(['Surgery', 'Non-surgery'], frameon=False, fontsize=8, loc='lower left')

# Panel B
plot_km(axes[0,1], df, 'surgery_vs_none',
        'B. CSS: Surgery vs Non-surgery', n_init, {'Surgery': C['seg'], 'Non-surgery': C['non_surg']})
axes[0,1].set_ylabel('Cancer-Specific Survival', fontsize=9)
axes[0,1].legend(['Surgery', 'Non-surgery'], frameon=False, fontsize=8, loc='lower left')

# Panel C — by surgery type
fits, labs = plot_km(axes[0,2], df, 'surgery_group', 'C. OS: By Surgery Type', {}, SURG_COLORS)
axes[0,2].legend(frameon=False, fontsize=6.5, loc='lower left')
add_risk_table(axes[0,2], fits, TIMES, labs, {lb: len(gr) for lb, gr in df.groupby('surgery_group') if len(gr)>=30})

# Panel D — by stage
fits_d, labs_d = [], []
for lb, gr in df.groupby('stage_label'):
    if len(gr) < 100: continue
    kmf = KaplanMeierFitter()
    kmf.fit(gr['surv_months'], gr['vital_dead'], label=lb)
    kmf.plot_survival_function(ax=axes[1,0], ci_show=False, lw=1.5)
    fits_d.append(kmf); labs_d.append(lb)
axes[1,0].set_title('D. OS: By AJCC Stage', fontweight='bold', fontsize=10)
axes[1,0].set_xlabel('Months', fontsize=9); axes[1,0].set_ylabel('Survival Probability', fontsize=9)
axes[1,0].set_xlim(0,120); axes[1,0].set_ylim(0,1.05)
axes[1,0].spines['top'].set_visible(False); axes[1,0].spines['right'].set_visible(False)
axes[1,0].legend(frameon=False, fontsize=8, loc='lower left')

# Panel E — surgery patients by age
for lb, gr in df[df['surgery_any']==1].groupby('age_band'):
    kmf = KaplanMeierFitter()
    kmf.fit(gr['surv_months'], gr['vital_dead'], label=lb)
    kmf.plot_survival_function(ax=axes[1,1], ci_show=False, lw=1.5)
axes[1,1].set_title('E. OS: Surgery Patients by Age', fontweight='bold', fontsize=10)
axes[1,1].set_xlabel('Months', fontsize=9); axes[1,1].set_ylabel('Survival Probability', fontsize=9)
axes[1,1].set_xlim(0,120); axes[1,1].set_ylim(0,1.05)
axes[1,1].spines['top'].set_visible(False); axes[1,1].spines['right'].set_visible(False)
axes[1,1].legend(frameon=False, fontsize=8, loc='lower left')

# Panel F — CSS by cancer type
for lb, gr in df.groupby('cancer_type'):
    if len(gr) < 500: continue
    kmf = KaplanMeierFitter()
    kmf.fit(gr['surv_months'], gr['css_dead'], label=lb)
    kmf.plot_survival_function(ax=axes[1,2], ci_show=False, lw=1.5)
axes[1,2].set_title('F. CSS: HCC vs ICC vs Other', fontweight='bold', fontsize=10)
axes[1,2].set_xlabel('Months', fontsize=9)
axes[1,2].set_ylabel('Cancer-Specific Survival', fontsize=9)
axes[1,2].set_xlim(0,120); axes[1,2].set_ylim(0,1.05)
axes[1,2].spines['top'].set_visible(False); axes[1,2].spines['right'].set_visible(False)
axes[1,2].legend(frameon=False, fontsize=8, loc='lower left')

plt.tight_layout(rect=[0, 0.03, 1, 1])
save(fig, 'Fig1_KM')
plt.close()

# ─── FIGURE 2: Forest Plot ─────────────────────────────────────────
print("Generating Fig2 — Forest plot...")
feats = ['age_c','male','married','race_nhb','race_nhapi','race_hispanic',
         'stage_2','stage_3','stage_4','grade_poor','is_icc',
         'surg_local_destruction','surg_segmental_resection',
         'surg_larger_resection','surg_transplant',
         'chemotherapy','radiation','cirrhosis']
labels_short = ['Age','Male','Married','Race: NHB','Race: NHAPI','Race: Hispanic',
                'Stage II','Stage III','Stage IV','Poor Grade','ICC Histology',
                'Local Destruction','Segmental Resection','Larger Resection','Transplant',
                'Chemotherapy','Radiation','Cirrhosis']

cox_df = df[['surv_months','css_dead']+feats].dropna()
m = CoxPHFitter(penalizer=0.01)
m.fit(cox_df[['surv_months','css_dead']+feats], 'surv_months', 'css_dead')

fig, ax = plt.subplots(figsize=(8, 7))
fig.patch.set_facecolor('white')
hrs, lo, hi = [], [], []
for c in feats:
    hr = np.exp(m.params_[c])
    cl, ch = get_ci(m, c)
    hrs.append(hr); lo.append(cl); hi.append(ch)

y = np.arange(len(feats))
# Color-code: surgery variables in blue, rest in gray
colors_forest = [C['seg'] if 'surg' in f or c == 'chemotherapy' or c == 'radiation' else '#333333' for f, c in zip(feats, feats)]

ax.errorbar(hrs, y, xerr=[np.array(hrs)-np.array(lo), np.array(hi)-np.array(hrs)],
            fmt='s', capsize=3, capthick=1, elinewidth=1, color='#333333',
            markerfacecolor=C['seg'], markeredgewidth=0)
ax.axvline(1, color='#c0392b', ls='--', lw=1, alpha=0.6)
ax.set_yticks(y)
ax.set_yticklabels(labels_short, fontsize=9)
ax.set_xlabel('Hazard Ratio (95% CI)', fontsize=10)
ax.set_xscale('log')
ax.set_xticks([0.1, 0.2, 0.5, 1, 2, 5])
ax.set_xlim(0.05, 8)
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
ax.axhline(y[-1]+0.5, color='#e0e0e0', lw=0.5)
ax.set_title('Figure 2. Multivariate Cox Regression — Cancer-Specific Survival', fontweight='bold', fontsize=10)

save(fig, 'Fig2_Forest')
plt.close()

# ─── FIGURE 3: HCC vs ICC ──────────────────────────────────────────
print("Generating Fig3 — HCC vs ICC comparison...")
fig, axes = plt.subplots(2, 3, figsize=(14, 11))
fig.patch.set_facecolor('white')

# A/B: Age-surgery HR by cancer type
for ct, ax_idx, color in [('HCC', (0,0), C['hcc']), ('ICC', (0,1), C['icc'])]:
    ax = axes[ax_idx]
    sub = df[df['cancer_type'] == ct]
    ages = [67, 72, 77, 85]
    for surg, label, ls in [('Segmental Resection', 'Segmental', '-'), ('Larger Resection', 'Larger', '--')]:
        hrs_row = []
        for lo, hi in [(64,70),(70,75),(75,80),(80,100)]:
            g = sub[(sub['age']>=lo)&(sub['age']<hi)]
            try:
                cx = CoxPHFitter()
                cx.fit(g[g['surgery_type'].isin(['None',surg])][['surv_months','css_dead','surgery_any']],
                       'surv_months','css_dead')
                h = np.exp(cx.params_['surgery_any'])
            except: h = np.nan
            hrs_row.append(h)
        ax.plot(ages, hrs_row, 'o-', lw=1.5, ls=ls, label=label, color=color, markersize=5)
    ax.axhline(1, color=C['ref_gray'], ls='-', lw=0.5)
    ax.set_xlabel('Age', fontsize=9); ax.set_ylabel('Surgery HR', fontsize=9)
    ax.set_title(f'{ct}: Segmental vs Larger Resection', fontweight='bold', fontsize=10)
    ax.set_ylim(0, 0.7); ax.legend(frameon=False, fontsize=8)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

# C: OS by surgery type HCC vs ICC
ax = axes[0,2]
kmf = KaplanMeierFitter()
for ct, color in [('HCC', C['hcc']), ('ICC', C['icc'])]:
    kmf.fit(df[df['cancer_type']==ct]['surv_months'], df[df['cancer_type']==ct]['vital_dead'], label=ct)
    kmf.plot_survival_function(ax=ax, ci_show=False, lw=1.5, color=color)
ax.set_title('C. Surgery Survival: HCC vs ICC', fontweight='bold', fontsize=10)
ax.set_xlabel('Months', fontsize=9); ax.set_ylabel('Survival Probability', fontsize=9)
ax.set_xlim(0, 60); ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
ax.legend(frameon=False, fontsize=8)

# D: Surgery distribution
ax = axes[1,0]
surg_order = ['None','Local Destruction','Segmental Resection','Larger Resection','Transplant']
surg_labels_short = ['None','Local','Segmental','Larger','Transplant']
x = np.arange(len(surg_order)); w = 0.35
for i, ct in enumerate(['HCC','ICC']):
    pcts = []
    for s in surg_order:
        n = len(df[(df['cancer_type']==ct)&(df['surgery_type']==s)])
        total = len(df[df['cancer_type']==ct])
        pcts.append(n/total*100)
    ax.bar(x + i*w - w/2, pcts, w, label=ct, color=C['hcc'] if ct=='HCC' else C['icc'], alpha=0.85)
ax.set_xticks(x); ax.set_xticklabels(surg_labels_short, fontsize=8)
ax.set_ylabel('%', fontsize=9); ax.set_title('D. Surgery Distribution', fontweight='bold', fontsize=10)
ax.legend(frameon=False, fontsize=8); ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

# E-F: Stage-specific HRs by cancer type
markers = {'Local_Destruction': ('o', '-', C['local']), 'Segmental_Resection': ('s', '-', C['seg']),
           'Larger_Resection': ('^', '--', C['larger']), 'Transplant': ('D', '-', C['txp'])}
for ci, ct in enumerate(['HCC','ICC']):
    ax = axes[1, ci+1]
    sub = df[df['cancer_type']==ct]
    for surg, (mkr, ls, clr) in markers.items():
        hrs_stage = []
        for st in [1,2,3,4]:
            g = sub[(sub['stage']==st)&(sub['surgery_type'].isin(['None',surg]))]
            try:
                cx = CoxPHFitter()
                cx.fit(g[['surv_months','css_dead','surgery_any']], 'surv_months','css_dead')
                h = np.exp(cx.params_['surgery_any'])
            except: h = np.nan
            hrs_stage.append(h)
        ax.plot([1,2,3,4], hrs_stage, marker=mkr, ls=ls, lw=1.5, markersize=7, label=surg.replace('_',' '), color=clr)
    ax.set_xlabel('AJCC Stage', fontsize=9); ax.set_ylabel('Surgery HR', fontsize=9)
    ax.set_title(f'{ct}: Stage-Specific HR', fontweight='bold', fontsize=10)
    ax.legend(frameon=False, fontsize=6.5); ax.axhline(1, color=C['ref_gray'], ls=':')
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

plt.tight_layout()
save(fig, 'Fig3_HCCvsICC')
plt.close()

print("\nDone — All JAMA-style figures saved to 04_Manuscript/figures/")
print("Files: Fig1_KM.png/pdf/tiff, Fig2_Forest.png/pdf/tiff, Fig3_HCCvsICC.png/pdf/tiff")
