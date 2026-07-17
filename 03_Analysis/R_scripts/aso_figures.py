"""ASO Figure Polish — Professional publication-ready styling
- Consistent ASO-spec: 174mm @ 300dpi TIFF LZW RGB
- Professional restrained palette (ASO surgical audience)
- Clean white backgrounds, no grid clutter
- Distinct colors for all 5 surgery groups
"""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import multivariate_logrank_test
from PIL import Image
import os, warnings; warnings.filterwarnings('ignore')
os.environ['PYTHONWARNINGS'] = 'ignore'

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.size': 8,
    'axes.titlesize': 9, 'axes.labelsize': 8,
    'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5,
    'legend.fontsize': 7, 'figure.dpi': 300,
    'axes.linewidth': 0.5, 'xtick.major.width': 0.4, 'ytick.major.width': 0.4,
})

BASE = r'D:\Researching\SEER\hepatobiliary cancer'
FIG = os.path.join(BASE, '04_Manuscript', 'figures')
os.chdir(BASE)
os.makedirs(FIG, exist_ok=True)

C = {
    'non_surg': '#7f8c8d', 'local': '#e67e22', 'seg': '#2980b9',
    'larger': '#8e44ad', 'txp': '#27ae60', 'hcc': '#2980b9',
    'icc': '#8e44ad', 'ref': '#bdc3c7', 'head': '#1a1a1a',
}
SURG = {'None':C['non_surg'], 'Local Destruction':C['local'],
        'Segmental Resection':C['seg'], 'Larger Resection':C['larger'],
        'Transplant':C['txp'], 'Other': C['ref']}
STAGE_C = {'I':'#3498db','II':'#2980b9','III':'#1f6dad','IV':'#1a5276'}
AGE_C  = {'65-69':'#e67e22','70-74':'#d35400','75-79':'#c0392b','80+':'#922b21'}

df = pd.read_csv(r'02_Data\cleaned\hepatobiliary_elderly_clean.csv')
df['surgery_type'] = df['surgery_type'].fillna('None').str.replace('_',' ').replace({'Surgery NOS':'Other','Other Surgery':'Other'})
df['surgery_any'] = (df['surgery_type'] != 'None').astype(int)
df['surgery_vs_none'] = df['surgery_any'].map({1:'Surgery',0:'Non-surgery'})
df['age_band'] = pd.cut(df['age'], bins=[64,70,75,80,100], labels=['65-69','70-74','75-79','80+'])
df['stage_label'] = df['stage'].map({1:'I',2:'II',3:'III',4:'IV'}).fillna('?')
df['surgery_group'] = df['surgery_type'].where(df['surgery_type'].isin(
    ['None','Local Destruction','Segmental Resection','Larger Resection','Transplant']), 'Other')

WIDTH_IN = 6.85
DPI = 300

def save(fig, name):
    fig.savefig(os.path.join(FIG, name)+'.png', dpi=DPI, bbox_inches='tight', facecolor='white')
    fig.savefig(os.path.join(FIG, name)+'.pdf', bbox_inches='tight', facecolor='white')
    im = Image.open(os.path.join(FIG, name)+'.png').convert('RGB')
    im.save(os.path.join(FIG, name)+'.tiff', 'TIFF', compression='tiff_lzw', dpi=(DPI,DPI))
    print(f'  ✓ {name}')

def style_ax(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(0.4)
    ax.spines['bottom'].set_linewidth(0.4)
    ax.tick_params(width=0.4)
    ax.set_facecolor('white')

def km_plot(ax, data, group, time_col, event_col, title, colors, max_t=120, risk_step=24,
            ylabel='Survival Probability'):
    groups = data[group].unique()
    groups = [g for g in groups if len(data[data[group]==g]) >= 30]
    groups.sort()
    dur_all, ev_all, grp_all = [], [], []
    medians = {}
    for g in groups:
        sub = data[data[group]==g]
        c = colors.get(g, '#333')
        kmf = KaplanMeierFitter()
        kmf.fit(sub[time_col], sub[event_col], label=g)
        m = kmf.median_survival_time_
        label = f'{g} ({int(m)}mo)' if np.isfinite(m) else g
        medians[g] = m
        kmf.plot_survival_function(ax=ax, ci_show=True, ci_alpha=0.12, lw=1.3, color=c, label=label)
        dur_all.extend(sub[time_col].tolist())
        ev_all.extend(sub[event_col].tolist())
        grp_all.extend([g]*len(sub))
    ax.axhline(0.5, color='#999', ls=':', lw=0.4, alpha=0.5, zorder=0)
    for g in groups:
        m = medians[g]
        if np.isfinite(m):
            c = colors.get(g, '#333')
            ax.plot([m, m], [0, 0.5], '--', lw=0.6, color=c, alpha=0.5, zorder=0)
    if len(groups) >= 2:
        lr = multivariate_logrank_test(dur_all, ev_all, grp_all)
        p = lr.p_value
        txt = f'Log-rank P {"<0.001" if p < 0.001 else f"={p:.3f}"}'
        ax.text(0.98, 0.98, txt, transform=ax.transAxes, fontsize=6.5,
                va='top', ha='right', color='#333', style='italic')
    ax.set_title(title, fontweight='bold', fontsize=8.5, color=C['head'], loc='left')
    ax.set_xlabel('Months', fontsize=7.5)
    ax.set_ylabel(ylabel, fontsize=7.5)
    ax.set_xlim(0, max_t); ax.set_ylim(0, 1.05)
    ax.set_xticks(np.arange(0, max_t+1, risk_step))
    style_ax(ax)
    ax.legend(frameon=False, fontsize=6.5, loc='lower left')

# ═══════════════════════════════════════════
# FIGURE 1: KM Curves — 2×3 panels
# ═══════════════════════════════════════════
print('Fig 1: Kaplan-Meier...')
fig, axes = plt.subplots(2, 3, figsize=(WIDTH_IN, 5.8))
fig.patch.set_facecolor('white')
plt.subplots_adjust(hspace=0.45, wspace=0.30, left=0.07, right=0.97, top=0.93, bottom=0.14)

# Bottom row: rotate x-labels
for i in range(3):
    axes[1,i].tick_params(axis='x', rotation=35)

km_plot(axes[0,0], df, 'surgery_vs_none', 'surv_months', 'vital_dead',
    'A. OS: Surgery vs No Surgery',
    {'Surgery':C['seg'], 'Non-surgery':C['non_surg']})

km_plot(axes[0,1], df, 'surgery_vs_none', 'surv_months', 'css_dead',
    'B. CSS: Surgery vs No Surgery',
    {'Surgery':C['seg'], 'Non-surgery':C['non_surg']}, ylabel='Cancer-Specific Survival')

km_plot(axes[0,2], df, 'surgery_group', 'surv_months', 'vital_dead',
    'C. OS by Surgery Type', SURG)

km_plot(axes[1,0], df, 'stage_label', 'surv_months', 'vital_dead',
    'D. OS by AJCC Stage', STAGE_C)

km_plot(axes[1,1], df[df['surgery_any']==1], 'age_band', 'surv_months', 'vital_dead',
    'E. OS by Age Group', AGE_C)

km_plot(axes[1,2], df, 'cancer_type', 'surv_months', 'css_dead',
    'F. HCC vs ICC', {'HCC':C['hcc'], 'ICC':C['icc']}, ylabel='Cancer-Specific Survival')

save(fig, 'Fig1_KM')
plt.close()

# ═══════════════════════════════════════════
# FIGURE 2: Forest Plot — JAMA style
# ═══════════════════════════════════════════
print('Fig 2: Forest plot...')
feats = ['age_c','male','married','race_nhb','race_nhapi','race_hispanic',
         'stage_2','stage_3','stage_4','grade_poor','is_icc',
         'surg_local_destruction','surg_segmental_resection',
         'surg_larger_resection','surg_transplant',
         'chemotherapy','radiation','cirrhosis']
labels_short = ['Age','Male','Married','Race: NHB','Race: NHAPI','Race: Hispanic',
                'Stage II','Stage III','Stage IV','Poor Grade','ICC',
                'Local Destruction','Segmental Resection','Larger Resection',
                'Transplant','Chemotherapy','Radiation','Cirrhosis']

cox_df = df[['surv_months','css_dead']+feats].dropna()
m = CoxPHFitter(penalizer=0.01)
m.fit(cox_df[['surv_months','css_dead']+feats], 'surv_months', 'css_dead')

hrs, lo, hi = [], [], []
for c in feats:
    hr = np.exp(m.params_[c])
    try:
        cl = np.exp(m.confidence_intervals_.loc[c,'95% lower-bound'])
        ch = np.exp(m.confidence_intervals_.loc[c,'95% upper-bound'])
    except:
        cl, ch = np.exp(m.params_[c]-1.96*m.standard_errors_[c]), np.exp(m.params_[c]+1.96*m.standard_errors_[c])
    hrs.append(hr); lo.append(cl); hi.append(ch)

y_var = np.arange(len(feats))

fig, ax = plt.subplots(figsize=(WIDTH_IN, 5.8))
fig.patch.set_facecolor('white')

cat_bands = [(0, 5, '#f5f6fa'), (6, 10, '#fafbfc'), (11, 14, '#f5f6fa'), (15, 17, '#fafbfc')]
for start, end, color in cat_bands:
    ax.axhspan(start-0.5, end+0.5, color=color, zorder=0)

colors_f = []
for i, c in enumerate(feats):
    if hrs[i] < 1:
        colors_f.append('#27ae60')  # green: protective
    else:
        colors_f.append('#c0392b')  # red: risk

ax.axvline(1, color='#c0392b', ls='--', lw=0.6, alpha=0.4, zorder=0)
for i in range(len(feats)):
    ax.plot([lo[i], hi[i]], [y_var[i], y_var[i]], '-', lw=2, color=colors_f[i], zorder=1)
    ax.plot(hrs[i], y_var[i], 'D', markersize=4.5, color=colors_f[i],
            markeredgecolor='white', markeredgewidth=0.3, zorder=2)

for i in range(len(feats)):
    ax.text(6, y_var[i], f'  {hrs[i]:.2f}  ({lo[i]:.2f}–{hi[i]:.2f})',
            va='center', fontsize=8, color=colors_f[i], fontfamily='monospace')

ax.set_yticks(y_var)
ax.set_yticklabels(labels_short, fontsize=8.5)
ax.set_xlabel('Hazard Ratio (95% CI)', fontsize=9)
ax.set_xscale('log')
ax.set_xticks([0.1, 0.2, 0.5, 1, 2, 5])
ax.set_xticklabels(['0.1','0.2','0.5','1','2','5'], fontsize=7.5)
ax.set_xlim(0.05, 6.8)
ax.set_ylim(-1.8, len(feats)+0.5)
ax.text(0.12, -2.2, '← Protective (HR<1)', fontsize=7.5, color='#27ae60', alpha=0.7)
ax.text(3, -2.2, 'Risk (HR>1) →', fontsize=7.5, color='#c0392b', alpha=0.7, ha='right')
ax.text(6, -2.2, 'HR (95% CI)', fontsize=8, color='#444', va='center', ha='center',
        fontweight='bold', fontfamily='monospace')
style_ax(ax)
ax.set_title('Figure 2. Multivariate Cox Regression: Cancer-Specific Survival',
             fontweight='bold', fontsize=9, color=C['head'], loc='left')
plt.tight_layout()
save(fig, 'Fig2_Forest')
plt.close()

# ═══════════════════════════════════════════
# FIGURE 3: HCC vs ICC Comparison
# ═══════════════════════════════════════════
print('Fig 3: HCC vs ICC...')
fig = plt.figure(figsize=(WIDTH_IN, 6.5))
fig.patch.set_facecolor('white')
gs2 = GridSpec(2, 3, hspace=0.45, wspace=0.35, left=0.08, right=0.97, bottom=0.15, top=0.95)
axes3 = [fig.add_subplot(gs2[i, j]) for i in range(2) for j in range(3)]

for ct, ax_idx, color in [('HCC',0,C['hcc']),('ICC',1,C['icc'])]:
    ax = axes3[ax_idx]; sub = df[df['cancer_type']==ct]
    ages = [67, 72, 77, 85]
    age_labels = ['65-69','70-74','75-79','80+']
    for surg, lb, ls in [('Segmental Resection','Segmental','-'),('Larger Resection','Larger','--')]:
        hrs_row = []
        for lo_a, hi_a in [(64,70),(70,75),(75,80),(80,100)]:
            g = sub[(sub['age']>=lo_a)&(sub['age']<hi_a)]
            try:
                cx = CoxPHFitter(penalizer=0.01)
                cx.fit(g[g['surgery_type'].isin(['None',surg])][['surv_months','css_dead','surgery_any']],
                       'surv_months','css_dead')
                h = np.exp(cx.params_['surgery_any'])
            except: h = np.nan
            hrs_row.append(h)
        ax.plot(ages, hrs_row, 'o-', lw=1.2, ls=ls, label=lb, color=color, markersize=4)
    ax.axhline(1, color=C['ref'], ls='-', lw=0.5)
    ax.set_xticks(ages); ax.set_xticklabels(age_labels, fontsize=6.5)
    ax.set_xlabel('Age Group', fontsize=7.5); ax.set_ylabel('Surgery HR (CSS)', fontsize=7.5)
    ax.set_title(f'{chr(65+ax_idx)}. {ct} Age HR', fontweight='bold', fontsize=8, color=C['head'], loc='left')
    ax.set_ylim(0,0.7); ax.legend(frameon=False, fontsize=6.5)
    style_ax(ax)

ax = axes3[2]
for ct, color in [('HCC',C['hcc']),('ICC',C['icc'])]:
    sub = df[df['cancer_type']==ct]
    kmf = KaplanMeierFitter()
    kmf.fit(sub['surv_months'], sub['css_dead'], label=ct)
    kmf.plot_survival_function(ax=ax, ci_show=True, ci_alpha=0.15, lw=1.3, color=color)
ax.set_title('C. HCC vs ICC', fontweight='bold', fontsize=8.5, color=C['head'], loc='left')
ax.set_xlabel('Months', fontsize=7.5); ax.set_ylabel('Cancer-Specific Survival', fontsize=7.5)
ax.set_xlim(0,60); ax.set_xticks(np.arange(0,61,12))
style_ax(ax); ax.legend(frameon=False, fontsize=7)

ax = axes3[3]
surg_order = ['None','Local Destruction','Segmental Resection','Larger Resection','Transplant']
x = np.arange(len(surg_order)); w = 0.3
labels_bar = ['None','Local','Segmental','Larger','Transplant']
for i, ct in enumerate(['HCC','ICC']):
    pcts = [len(df[(df['cancer_type']==ct)&(df['surgery_type']==s)])/len(df[df['cancer_type']==ct])*100 for s in surg_order]
    ax.bar(x + i*w - w/2, pcts, w, label=ct, color=C['hcc'] if ct=='HCC' else C['icc'], alpha=0.85)
ax.set_xticks(x); ax.set_xticklabels(labels_bar, fontsize=6.5, rotation=45)
ax.set_ylabel('%', fontsize=7.5); ax.set_title('D. Surgery by Type', fontweight='bold', fontsize=8.5, color=C['head'], loc='left')
ax.legend(frameon=False, fontsize=7); style_ax(ax)

markers = {'Local_Destruction':('o','-',C['local']), 'Segmental_Resection':('s','-',C['seg']),
           'Larger_Resection':('^','--',C['larger']), 'Transplant':('D','-',C['txp'])}
for ci, ct in enumerate(['HCC','ICC']):
    ax = axes3[4+ci]; sub = df[df['cancer_type']==ct]
    for surg, (mkr, ls, clr) in markers.items():
        sname = surg.replace('_',' ')
        hrs_st = []
        for st in [1,2,3,4]:
            g = sub[(sub['stage']==st)&(sub['surgery_type'].isin(['None',sname]))]
            try:
                cx = CoxPHFitter(penalizer=0.01)
                cx.fit(g[['surv_months','css_dead','surgery_any']],'surv_months','css_dead')
                h = np.exp(cx.params_['surgery_any'])
            except: h = np.nan
            hrs_st.append(h)
        ax.plot([1,2,3,4], hrs_st, marker=mkr, ls=ls, lw=1.2, markersize=5,
                label=surg.replace('_',' ').replace('Local Destruction','Local Destr.').replace('Resection','Resec.'), color=clr)
    ax.set_xlabel('AJCC Stage', fontsize=7.5); ax.set_ylabel('Surgery HR', fontsize=7.5)
    ax.set_title(f'{chr(69+ci)}. {ct} Stage HR', fontweight='bold', fontsize=8, color=C['head'], loc='left')
    ax.legend(frameon=False, fontsize=5, loc='upper right'); ax.axhline(1, color=C['ref'], ls=':')
    ax.set_xticks([1,2,3,4])
    style_ax(ax)

save(fig, 'Fig3_HCCvsICC')
plt.close()

print()
print('─'*60)
print('ASO figures generated:')
print('  • Fig1: KM curves with median survival + log-rank P')
print('  • Fig2: Forest plot with category bands + HR column (JAMA style)')
print('  • Fig20: HCC vs ICC sub-analyses')
print('─'*60)
