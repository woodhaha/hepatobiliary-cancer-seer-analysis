"""Phases 3-5: Baseline Table + KM + Cox + Stratification — Hepatobiliary Cancer"""
import pandas as pd
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test
from scipy.stats import chi2_contingency, ttest_ind
import os, warnings
warnings.filterwarnings('ignore')

os.chdir(r'D:\Researching\SEER\hepatobiliary cancer')
os.makedirs('03_Analysis/figures', exist_ok=True)
os.makedirs('03_Analysis/outputs', exist_ok=True)

df = pd.read_csv(r'02_Data\cleaned\hepatobiliary_elderly_clean.csv')
print(f"Loaded: {len(df)} elderly")

# Fix surgery_type NaN
df['surgery_type'] = df['surgery_type'].fillna('None')
df['surgery_any'] = (df['surgery_type'] != 'None').astype(int)
df['surgery_vs_none'] = df['surgery_type'].apply(lambda x: 'Surgery' if x != 'None' else 'Non-surgery')
df['surgery_group'] = df['surgery_type'].map({
    'None': 'Non-surgery', 'Local_Destruction': 'Local Destruction',
    'Segmental_Resection': 'Segmental Resection', 'Larger_Resection': 'Larger Resection',
    'Transplant': 'Transplant', 'Other_Surgery': 'Other', 'Surgery_NOS': 'Other'})
df['stage_label'] = df['stage'].map({1:'Stage I',2:'Stage II',3:'Stage III',4:'Stage IV'})
df['age_band'] = pd.cut(df['age'], bins=[64,70,75,80,100], labels=['65-69','70-74','75-79','80+'])

# Helper for safe CI extraction
def get_ci(cph, col):
    ci = cph.confidence_intervals_
    cols = list(ci.columns)
    try:
        return np.exp(ci.loc[col, cols[0]]), np.exp(ci.loc[col, cols[1]])
    except:
        return np.nan, np.nan

ns = df[df['surgery_type']=='None']
sg = df[df['surgery_type']!='None']

with open('03_Analysis/outputs/02_survival_report.md', 'w', encoding='utf-8') as rpt:
    def p(*a, **kw):
        print(*a, **kw)
        print(*a, **kw, file=rpt)

    p("# Hepatobiliary Cancer Survival Analysis\n")
    p(f"**Cohort**: {len(df)} elderly (≥65), {len(ns)} Non-surgery, {len(sg)} Surgery\n")

    # ============================================================
    # TABLE 1
    # ============================================================
    p("## Table 1: Baseline Characteristics\n")
    p("| Variable | Non-surgery (N={}) | Surgery (N={}) | P-value |".format(len(ns), len(sg)))
    p("|----------|-------------------|----------------|---------|")

    baseline = [
        ('Age (mean±SD)', 'age', 'cont'),
        ('Male (%)', 'male', 'bin'),
        ('Married (%)', 'married', 'bin'),
        ('Stage I (%)', lambda d: (d['stage']==1).mean(), 'func'),
        ('Stage II (%)', lambda d: (d['stage']==2).mean(), 'func'),
        ('Stage III (%)', lambda d: (d['stage']==3).mean(), 'func'),
        ('Stage IV (%)', lambda d: (d['stage']==4).mean(), 'func'),
        ('Poor Grade (%)', 'grade_poor', 'bin'),
        ('ICC (%)', 'is_icc', 'bin'),
        ('Chemotherapy (%)', 'chemotherapy', 'bin'),
        ('Radiation (%)', 'radiation', 'bin'),
        ('Cirrhosis (%)', 'cirrhosis', 'bin'),
    ]

    for label, var, vtype in baseline:
        if vtype == 'cont':
            ns_v = f"{ns[var].mean():.1f}±{ns[var].std():.1f}"
            sg_v = f"{sg[var].mean():.1f}±{sg[var].std():.1f}"
            try: _, pv = ttest_ind(ns[var].dropna(), sg[var].dropna())
            except: pv = 1.0
        elif vtype == 'bin':
            ns_v = f"{ns[var].mean()*100:.1f}%"
            sg_v = f"{sg[var].mean()*100:.1f}%"
            try:
                tbl = pd.crosstab(df['surgery_type']=='None', df[var])
                _, pv, _, _ = chi2_contingency(tbl)
            except: pv = 1.0
        elif vtype == 'func':
            ns_v = f"{var(ns)*100:.1f}%"
            sg_v = f"{var(sg)*100:.1f}%"
            pv = 0.5
        ps = f"{pv:.3f}" if pv >= 0.001 else "<0.001"
        p(f"| {label} | {ns_v} | {sg_v} | {ps} |")

    p()

    # ============================================================
    # FIGURE 1: KM Curves
    # ============================================================
    p("## Figure 1: Kaplan-Meier\n")
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle('Hepatobiliary Cancer — Elderly (≥65)', fontsize=16, fontweight='bold')
    kmf = KaplanMeierFitter()

    # A: OS Surgery vs Non
    ax = axes[0,0]
    for lb, gr in df.groupby('surgery_vs_none'):
        kmf.fit(gr['surv_months'], gr['vital_dead'], label=lb)
        kmf.plot_survival_function(ax=ax, ci_show=False, lw=2)
    ax.set_title('A. OS: Surgery vs Non-surgery', fontweight='bold')
    ax.set(xlabel='Months', ylabel='Overall Survival', xlim=(0,120)); ax.legend(frameon=False)

    # B: CSS Surgery vs Non
    ax = axes[0,1]
    for lb, gr in df.groupby('surgery_vs_none'):
        kmf.fit(gr['surv_months'], gr['css_dead'], label=lb)
        kmf.plot_survival_function(ax=ax, ci_show=False, lw=2)
    ax.set_title('B. CSS: Surgery vs Non-surgery', fontweight='bold')
    ax.set(xlabel='Months', ylabel='Cancer-Specific Survival', xlim=(0,120)); ax.legend(frameon=False)

    # C: OS by surgery type
    ax = axes[0,2]
    colors = {'Non-surgery':'#e74c3c','Local Destruction':'#f39c12','Segmental Resection':'#2ecc71',
              'Larger Resection':'#3498db','Transplant':'#9b59b6'}
    for lb, gr in df.groupby('surgery_group'):
        if len(gr)>=30:
            kmf.fit(gr['surv_months'], gr['vital_dead'], label=lb)
            kmf.plot_survival_function(ax=ax, ci_show=False, lw=2, color=colors.get(lb,'gray'))
    ax.set_title('C. OS: By Surgery Type', fontweight='bold')
    ax.set(xlabel='Months', ylabel='Overall Survival', xlim=(0,120)); ax.legend(frameon=False, fontsize=8)

    # D: OS by stage
    ax = axes[1,0]
    for lb, gr in df.groupby('stage_label'):
        kmf.fit(gr['surv_months'], gr['vital_dead'], label=lb)
        kmf.plot_survival_function(ax=ax, ci_show=False, lw=2)
    ax.set_title('D. OS: By AJCC Stage', fontweight='bold')
    ax.set(xlabel='Months', ylabel='Overall Survival', xlim=(0,120)); ax.legend(frameon=False)

    # E: OS surgery patients by age
    ax = axes[1,1]
    for lb, gr in df[df['surgery_any']==1].groupby('age_band'):
        kmf.fit(gr['surv_months'], gr['vital_dead'], label=lb)
        kmf.plot_survival_function(ax=ax, ci_show=False, lw=2)
    ax.set_title('E. OS: Surgery Patients by Age', fontweight='bold')
    ax.set(xlabel='Months', ylabel='Overall Survival', xlim=(0,120)); ax.legend(frameon=False)

    # F: CSS by cancer type
    ax = axes[1,2]
    for lb, gr in df.groupby('cancer_type'):
        if len(gr)>=500:
            kmf.fit(gr['surv_months'], gr['css_dead'], label=lb)
            kmf.plot_survival_function(ax=ax, ci_show=False, lw=2)
    ax.set_title('F. CSS: HCC vs ICC vs Other', fontweight='bold')
    ax.set(xlabel='Months', ylabel='Cancer-Specific Survival', xlim=(0,120)); ax.legend(frameon=False)

    plt.tight_layout()
    fig.savefig('03_Analysis/figures/Fig1_KM.png', dpi=300, bbox_inches='tight')
    fig.savefig('03_Analysis/figures/Fig1_KM.pdf', bbox_inches='tight')
    plt.close()
    p("✓ Fig1 saved\n")

    # ============================================================
    # LOG-RANK
    # ============================================================
    p("## Log-Rank Tests\n")
    for oc, nm in [('vital_dead','OS'),('css_dead','CSS')]:
        r = logrank_test(ns['surv_months'], sg['surv_months'], ns[oc], sg[oc])
        p(f"{nm} Surgery vs Non: χ²={r.test_statistic:.1f}, p={r.p_value:.4f}")
    p()

    # ============================================================
    # COX REGRESSION
    # ============================================================
    p("## Cox Regression\n")
    feats = ['age_c','male','married','race_nhb','race_nhapi','race_hispanic',
             'stage_2','stage_3','stage_4','grade_poor','is_icc',
             'surg_local_destruction','surg_segmental_resection',
             'surg_larger_resection','surg_transplant',
             'chemotherapy','radiation','cirrhosis']
    cox_df = df[['surv_months','css_dead']+feats].dropna()

    # Univariate
    p("### Univariate CSS Cox\n| Variable | HR | 95% CI | P |")
    p("|---|---|---|---|")
    uva = []
    for c in feats:
        try:
            m = CoxPHFitter(); m.fit(cox_df[['surv_months','css_dead',c]],'surv_months','css_dead')
            hr = np.exp(m.params_[c]); ci_l, ci_h = get_ci(m, c); pv = m.summary.loc[c,'p']
            uva.append({'var':c,'HR':hr,'CI_low':ci_l,'CI_high':ci_h,'p':pv})
            p(f"| {c} | {hr:.2f} | [{ci_l:.2f}-{ci_h:.2f}] | {pv:.3f} |")
        except Exception as e: pass
    p()

    # Multivariate
    p("### Multivariate CSS Cox\n| Variable | HR | 95% CI | P |")
    p("|---|---|---|---|")
    try:
        m = CoxPHFitter(penalizer=0.01)
        m.fit(cox_df[['surv_months','css_dead']+feats],'surv_months','css_dead')
        for c in feats:
            if c in m.params_.index:
                hr = np.exp(m.params_[c]); ci_l, ci_h = get_ci(m, c); pv = m.summary.loc[c,'p']
                ps = f"<0.001" if pv < 0.001 else f"{pv:.3f}"
                p(f"| {c} | {hr:.2f} | [{ci_l:.2f}-{ci_h:.2f}] | {ps} |")
        p(f"\nC-index={m.concordance_index_:.3f}")
    except Exception as e:
        p(f"Cox failed: {e}")
    p()

    # ============================================================
    # STRATIFICATION (Paper Table 3 style)
    # ============================================================
    p("## Stratification Analysis\n")
    p("| Strata | NoSurg N | LD HR[CI] | SR HR[CI] | LR HR[CI] | LT HR[CI] |")
    p("|--------|---------|-----------|-----------|-----------|-----------|")

    for strata_var in ['age_band','stage_label','cancer_type']:
        p(f"**{strata_var}** | | | | | |")
        for st, gr in df.groupby(strata_var):
            if len(gr)<100: continue
            ref_n = (gr['surgery_type']=='None').sum()
            row = [str(st), str(ref_n)]
            for s in ['Local_Destruction','Segmental_Resection','Larger_Resection','Transplant']:
                sub = gr[gr['surgery_type'].isin(['None',s])]
                try:
                    cx = CoxPHFitter()
                    cx.fit(sub[['surv_months','css_dead','surgery_any']],'surv_months','css_dead')
                    hr = np.exp(cx.params_['surgery_any'])
                    ci_l, ci_h = get_ci(cx, 'surgery_any')
                    row.append(f"{hr:.2f}[{ci_l:.2f}-{ci_h:.2f}]")
                except: row.append("—")
            p("| "+" | ".join(row)+" |")

    p()

    # ============================================================
    # KEY INSIGHTS
    # ============================================================
    p("## Key Insights\n")
    p("### Median Survival by Surgery Type")
    for st, gr in df.groupby('surgery_group'):
        p(f"  {st}: median OS={gr['surv_months'].median():.0f}m, mortality={gr['vital_dead'].mean()*100:.1f}%")

    p("\n### Age-Surgery Benefit")
    for ab, gr in df.groupby('age_band'):
        ns_d = gr[gr['surgery_type']=='None']['vital_dead'].mean()
        s_d = gr[gr['surgery_type']!='None']['vital_dead'].mean()
        p(f"  {ab}: mortality gap = {(ns_d-s_d)*100:.1f}% (NonSurg {ns_d*100:.0f}% vs Surg {s_d*100:.0f}%)")

    p("\n### HCC vs ICC")
    for ct, gr in df.groupby('cancer_type'):
        p(f"  {ct}: median OS={gr['surv_months'].median():.0f}m, surgery rate={gr['surgery_any'].mean()*100:.1f}%")

    p("\n### Surgery Benefit by Stage (CSS HR)")
    for st, gr in df.groupby('stage'):
        if len(gr)<200: continue
        try:
            cx = CoxPHFitter()
            cx.fit(gr[['surv_months','css_dead','surgery_any']],'surv_months','css_dead')
            hr = np.exp(cx.params_['surgery_any'])
            p(f"  Stage {int(st)}: surgery CSS HR={hr:.2f}, N={len(gr)}")
        except: pass

    # ============================================================
    # FIGURE 2: FOREST PLOT
    # ============================================================
    p("\n## Figure 2: Forest Plot\n")
    try:
        fig, ax = plt.subplots(figsize=(10, 8))
        hrs, lo, hi = [], [], []
        for c in feats:
            if c in m.params_.index:
                hrs.append(np.exp(m.params_[c]))
                cl, ch = get_ci(m, c)
                lo.append(cl); hi.append(ch)
            else:
                hrs.append(np.nan); lo.append(np.nan); hi.append(np.nan)

        y = range(len(feats))
        ax.errorbar(hrs, y, xerr=[np.array(hrs)-np.array(lo), np.array(hi)-np.array(hrs)],
                    fmt='o', capsize=3, color='steelblue')
        ax.axvline(1, color='red', ls='--', alpha=0.5)
        ax.set_yticks(y); ax.set_yticklabels(feats, fontsize=9)
        ax.set_xlabel('Hazard Ratio (95% CI)'); ax.set_xscale('log')
        ax.set_title('Multivariate CSS Cox — Elderly Hepatobiliary Cancer', fontweight='bold')
        fig.savefig('03_Analysis/figures/Fig2_Forest.png', dpi=300, bbox_inches='tight')
        fig.savefig('03_Analysis/figures/Fig2_Forest.pdf', bbox_inches='tight')
        plt.close()
        p("✓ Fig2 saved")
    except Exception as e:
        p(f"Forest plot failed: {e}")

print("\nDONE — Report: 03_Analysis/outputs/02_survival_report.md")
