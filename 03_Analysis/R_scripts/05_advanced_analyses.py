"""Phase 9: Advanced Analyses — SHAP, Competing Risk, RMST, Age-Surgery Spline, Nomogram"""
import pandas as pd
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=['#0072B2','#E69F00','#009E73','#CC79A7','#56B4E9','#F0E442','#000000'])
from PIL import Image
import os, warnings, json
warnings.filterwarnings('ignore')

os.chdir(r'D:\Researching\SEER\hepatobiliary cancer')
os.makedirs('03_Analysis/figures', exist_ok=True)
os.makedirs('03_Analysis/outputs', exist_ok=True)

FIGDIR = '04_Manuscript/figures'
W = 6.85
DPI = 300
os.makedirs(FIGDIR, exist_ok=True)

def _sax(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(0.35)
    ax.spines['bottom'].set_linewidth(0.35)
    ax.tick_params(width=0.35)
    ax.set_facecolor('white')

_ASO_RC = {
    'font.family': 'sans-serif', 'font.size': 7,
    'axes.titlesize': 7.5, 'axes.labelsize': 7,
    'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5,
    'legend.fontsize': 6, 'figure.dpi': 300,
    'axes.linewidth': 0.4, 'xtick.major.width': 0.35, 'ytick.major.width': 0.35,
}

df = pd.read_csv(r'02_Data\cleaned\hepatobiliary_elderly_clean.csv')
df['surgery_type'] = df['surgery_type'].fillna('None')
df['surgery_any'] = (df['surgery_type'] != 'None').astype(int)
df['age_band'] = pd.cut(df['age'], bins=[64,70,75,80,100], labels=['65-69','70-74','75-79','80+'])
df['stage_label'] = df['stage'].map({1:'Stage I',2:'Stage II',3:'Stage III',4:'Stage IV'})
df['tumor_size'] = df['tumor_size'].fillna(df['tumor_size'].median())
df['year'] = df.get('year', 2010)  # fallback

features = ['age_c','male','married','race_nhb','race_nhapi','race_hispanic',
            'stage_2','stage_3','stage_4','grade_poor','is_icc',
            'chemotherapy','radiation','cirrhosis','income_10k','tumor_size','surgery_any']

with open('03_Analysis/outputs/04_advanced_report.md', 'w', encoding='utf-8') as rpt:
    def p(*a, **kw):
        print(*a, **kw)
        print(*a, **kw, file=rpt)

    p("# Advanced Analyses — Hepatobiliary Cancer\n")

    # ============================================================
    # 1. SHAP ANALYSIS
    # ============================================================
    p("## 1. SHAP Feature Analysis\n")
    try:
        from sksurv.ensemble import RandomSurvivalForest
        from sksurv.util import Surv
        import shap

        ml_df = df[['surv_months','vital_dead']+features].dropna().copy()
        X = ml_df[features].values.astype(np.float32)
        y = Surv.from_arrays(event=ml_df['vital_dead'].values.astype(bool),
                             time=ml_df['surv_months'].values.astype(np.float64))

        # Train RSF
        rsf = RandomSurvivalForest(n_estimators=100, min_samples_split=50,
                                   min_samples_leaf=20, max_features='sqrt',
                                   n_jobs=-1, random_state=42)
        rsf.fit(X[:30000], y[:30000])

        # SHAP on XGBoost-style
        import xgboost as xgb
        dtrain = xgb.DMatrix(X[:30000])
        dtrain.set_float_info('label', y['time'][:30000])
        dtrain.set_float_info('label_lower_bound', y['time'][:30000])
        dtrain.set_float_info('label_upper_bound',
            np.where(y['event'][:30000], y['time'][:30000], np.inf))
        params = {'objective':'survival:cox','eval_metric':'cox-nloglik',
                  'max_depth':4,'learning_rate':0.05,'subsample':0.8,
                  'colsample_bytree':0.8,'tree_method':'hist','seed':42}
        xgb_model = xgb.train(params, dtrain, num_boost_round=150)

        # SHAP TreeExplainer
        explainer = shap.TreeExplainer(xgb_model)
        shap_values = explainer.shap_values(X[:1000])

        plt.rcParams.update(**_ASO_RC)
        fig, ax = plt.subplots(figsize=(10, 8))
        shap_feat_labels = {
            'age_c':'Age','male':'Male','married':'Married',
            'race_nhb':'Race: NHB','race_nhapi':'Race: NHAPI','race_hispanic':'Race: Hispanic',
            'stage_2':'Stage II','stage_3':'Stage III','stage_4':'Stage IV',
            'grade_poor':'Poor Grade','is_icc':'ICC',
            'chemotherapy':'Chemotherapy','radiation':'Radiation',
            'cirrhosis':'Cirrhosis','income_10k':'Income','tumor_size':'Tumor Size','surgery_any':'Surgery',
        }
        shap_names = [shap_feat_labels.get(f, f) for f in features]
        shap.summary_plot(shap_values, X[:1000], feature_names=shap_names,
                         show=False, max_display=15)
        _sax(ax)
        fig.savefig('03_Analysis/figures/FigS1_SHAP.png', dpi=300, bbox_inches='tight')
        fig.set_size_inches(W, fig.get_size_inches()[1] * W / fig.get_size_inches()[0])
        fig.savefig(os.path.join(FIGDIR, 'FigS1_SHAP.png'), dpi=DPI, bbox_inches='tight', facecolor='white')
        fig.savefig(os.path.join(FIGDIR, 'FigS1_SHAP.pdf'), bbox_inches='tight', facecolor='white')
        Image.open(os.path.join(FIGDIR, 'FigS1_SHAP.png')).convert('RGB').save(os.path.join(FIGDIR, 'FigS1_SHAP.tiff'), 'TIFF', compression='tiff_lzw', dpi=(DPI,DPI))
        plt.close()
        p("✓ SHAP beeswarm saved")
    except Exception as e:
        p(f"SHAP failed: {e}")

    # ============================================================
    # 2. COMPETING RISK (CSS vs Other-Cause Death)
    # ============================================================
    p("\n## 2. Competing Risk Analysis\n")

    # Status codes: 0=alive, 1=CSS death, 2=other death
    df['comp_status'] = 0
    df.loc[df['css_dead']==1, 'comp_status'] = 1
    df.loc[(df['vital_dead']==1) & (df['css_dead']==0), 'comp_status'] = 2

    p("| Cause | N | % of Deaths |")
    p("|---|---|---|")
    p(f"| Cancer-Specific Death | {df['css_dead'].sum()} | {df['css_dead'].sum()/df['vital_dead'].sum()*100:.1f}% |")
    p(f"| Other-Cause Death | {(df['vital_dead']==1).sum() - df['css_dead'].sum()} | {(1-df['css_dead'].sum()/df['vital_dead'].sum())*100:.1f}% |")
    p(f"| Alive | {(df['vital_dead']==0).sum()} | — |")

    # CIF by surgery group — manual Nelson-Aalen style
    from lifelines import AalenJohansenFitter

    plt.rcParams.update(**_ASO_RC)
    fig, axes = plt.subplots(1, 3, figsize=(W, 3.5))

    # By surgery
    for i, (groupby, title) in enumerate([
        ('surgery_any', 'A. Surgery vs Non-surgery'),
        ('age_band', 'B. By Age'),
        ('stage_label', 'C. By Stage')
    ]):
        ax = axes[i]
        for name, grp in df.groupby(groupby):
            if len(grp) < 100: continue
            ajf = AalenJohansenFitter()
            lbl_map = {1:'Surgery',0:'Non-surgery'} if groupby=='surgery_any' else {}
            label = lbl_map.get(name, str(name))
            ajf.fit(grp['surv_months'], grp['comp_status'],
                   event_of_interest=1, label=label)
            ajf.plot(ax=ax, ci_show=False, lw=1.2)
        ax.set_title(title, fontweight='bold', fontsize=7.5, loc='left')
        ax.set_ylabel('CSS Cumulative Incidence', fontsize=7)
        ax.set_xlabel('Months', fontsize=7)
        ax.set_xlim(0, 120)
        ax.set_xticks(np.arange(0, 121, 24))
        ax.legend(frameon=False, fontsize=5.5)
        _sax(ax)

    fig.savefig('03_Analysis/figures/FigS2_CompetingRisk.png', dpi=300, bbox_inches='tight')
    fig.savefig('03_Analysis/figures/FigS2_CompetingRisk.pdf', bbox_inches='tight')
    fig.set_size_inches(W, 3.5)
    fig.savefig(os.path.join(FIGDIR, 'FigS2_CompetingRisk.png'), dpi=DPI, bbox_inches='tight', facecolor='white')
    fig.savefig(os.path.join(FIGDIR, 'FigS2_CompetingRisk.pdf'), bbox_inches='tight', facecolor='white')
    Image.open(os.path.join(FIGDIR, 'FigS2_CompetingRisk.png')).convert('RGB').save(os.path.join(FIGDIR, 'FigS2_CompetingRisk.tiff'), 'TIFF', compression='tiff_lzw', dpi=(DPI,DPI))
    plt.close()
    p("✓ Competing risk CIF saved")

    # ============================================================
    # 3. RMST (Restricted Mean Survival Time)
    # ============================================================
    p("\n## 3. RMST at 60 Months\n")
    from lifelines.utils import restricted_mean_survival_time

    p("| Group | RMST (months) | SE |")
    p("|---|---|---|")

    for groupby, label in [('surgery_type','Surgery Type'), ('stage_label','Stage'), ('age_band','Age')]:
        p(f"**{label}** | | |")
        for name, grp in df.groupby(groupby):
            try:
                rmst = restricted_mean_survival_time(
                    grp['surv_months'].values,
                    grp['vital_dead'].values.astype(bool),
                    t=60)
                p(f"| {name} | {rmst:.1f} | — |")
            except: pass
        p()

    # ============================================================
    # 4. AGE-SURGERY BENEFIT SPLINE
    # ============================================================
    p("## 4. Age-Surgery Benefit Curve\n")
    from lifelines import CoxPHFitter

    # Compute surgery HR in 2-year age bins
    ages = np.arange(65, 91, 2)
    hrs, ci_lo, ci_hi = [], [], []

    for a_min in ages:
        sub = df[(df['age'] >= a_min) & (df['age'] < a_min + 2)]
        if len(sub) < 200: continue
        try:
            cph = CoxPHFitter()
            cph.fit(sub[['surv_months','css_dead','surgery_any','stage_2','stage_3','stage_4',
                         'grade_poor','is_icc','male','chemotherapy']],
                    'surv_months','css_dead')
            hr = np.exp(cph.params_['surgery_any'])
            ci = cph.confidence_intervals_.loc['surgery_any']
            ci_l, ci_h = np.exp(ci.iloc[0]), np.exp(ci.iloc[1])
            hrs.append(hr); ci_lo.append(ci_l); ci_hi.append(ci_h)
        except:
            hrs.append(np.nan); ci_lo.append(np.nan); ci_hi.append(np.nan)

    plt.rcParams.update(**_ASO_RC)
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(ages[:len(hrs)], hrs, 'o-', color='steelblue', lw=2)
    ax.fill_between(ages[:len(hrs)], ci_lo, ci_hi, alpha=0.2, color='steelblue')
    ax.axhline(y=1, color='red', ls='--', alpha=0.5)

    # Compare with paper: segmental vs larger
    for surg_type in ['Segmental_Resection','Larger_Resection']:
        hrs2 = []
        for a_min in ages:
            sub = df[(df['age'] >= a_min) & (df['age'] < a_min + 2)]
            sub = sub[sub['surgery_type'].isin(['None', surg_type])]
            if len(sub) < 100: continue
            try:
                cph = CoxPHFitter()
                cph.fit(sub[['surv_months','css_dead','surgery_any','stage_2','stage_3','stage_4',
                             'grade_poor','is_icc','male']],
                        'surv_months','css_dead')
                hrs2.append(np.exp(cph.params_['surgery_any']))
            except: hrs2.append(np.nan)

        label = 'Segmental' if 'Segmental' in surg_type else 'Larger Resection'
        ax.plot(ages[:len(hrs2)], hrs2, 's--', lw=1.5,
               label=f'{label} vs None', alpha=0.7)

    ax.set_xlabel('Age'); ax.set_ylabel('Surgery CSS HR (adjusted)')
    ax.set_title('Age-Dependent Surgery Benefit in Hepatobiliary Cancer', fontweight='bold')
    ax.legend(frameon=False)
    ax.set_ylim(0, 0.8)
    _sax(ax)

    fig.savefig('03_Analysis/figures/FigS3_AgeSurgeryBenefit.png', dpi=300, bbox_inches='tight')
    fig.savefig('03_Analysis/figures/FigS3_AgeSurgeryBenefit.pdf', bbox_inches='tight')
    fig.set_size_inches(W, fig.get_size_inches()[1] * W / fig.get_size_inches()[0])
    fig.savefig(os.path.join(FIGDIR, 'FigS3_AgeSurgeryBenefit.png'), dpi=DPI, bbox_inches='tight', facecolor='white')
    fig.savefig(os.path.join(FIGDIR, 'FigS3_AgeSurgeryBenefit.pdf'), bbox_inches='tight', facecolor='white')
    Image.open(os.path.join(FIGDIR, 'FigS3_AgeSurgeryBenefit.png')).convert('RGB').save(os.path.join(FIGDIR, 'FigS3_AgeSurgeryBenefit.tiff'), 'TIFF', compression='tiff_lzw', dpi=(DPI,DPI))
    plt.close()
    p("✓ Age-surgery benefit spline saved")

    # ============================================================
    # 5. HEPATOBILIARY-SPECIFIC INSIGHTS
    # ============================================================
    p("\n## 5. Novel Insights — Hepatobiliary Cancer\n")

    # Insight A: HCC vs ICC surgery patterns
    p("### A. HCC vs ICC Treatment Gap\n")
    for ct in ['HCC','ICC']:
        sub = df[df['cancer_type']==ct]
        p(f"**{ct}**: N={len(sub)}, Surgery={sub['surgery_any'].mean()*100:.1f}%")
        for st in ['None','Local_Destruction','Segmental_Resection','Larger_Resection','Transplant']:
            n = (sub['surgery_type']==st).sum()
            med = sub[sub['surgery_type']==st]['surv_months'].median()
            p(f"  {st}: N={n}, median OS={med:.0f}m")

    # Insight B: Cirrhosis paradox
    p("\n### B. Cirrhosis Paradox\n")
    for cirr in [0, 1]:
        sub = df[df['cirrhosis']==cirr]
        p(f"Cirrhosis={cirr}: surgery rate={sub['surgery_any'].mean()*100:.1f}%, "
          f"median OS={sub['surv_months'].median():.0f}m, mortality={sub['vital_dead'].mean()*100:.1f}%")

    # Insight C: Yearly trends
    p("\n### C. Temporal Trends (2000-2022)\n")
    p("| Year | N | Surgery % | Transplant % | Median OS |")
    p("|---|---|---|---|---|")
    for yr in range(2000, 2023):
        sub = df[df['year']==yr]
        if len(sub) < 100: continue
        p(f"| {yr} | {len(sub)} | {sub['surgery_any'].mean()*100:.0f}% | "
          f"{sub['surgery_type'].eq('Transplant').mean()*100:.1f}% | "
          f"{sub['surv_months'].median():.0f}m |")

    # Trend figure
    yrs = range(2000, 2023)
    trends = {'Surgery Rate': [], 'Transplant Rate': [], 'Median OS': []}
    for yr in yrs:
        sub = df[df['year']==yr]
        trends['Surgery Rate'].append(sub['surgery_any'].mean()*100 if len(sub)>100 else np.nan)
        trends['Transplant Rate'].append(sub['surgery_type'].eq('Transplant').mean()*100 if len(sub)>100 else np.nan)
        trends['Median OS'].append(sub['surv_months'].median() if len(sub)>100 else np.nan)

    plt.rcParams.update(**_ASO_RC)
    fig, ax1 = plt.subplots(figsize=(10, 5))
    ax1.plot(yrs, trends['Surgery Rate'], 'o-', color='steelblue', label='Surgery Rate %')
    ax1.plot(yrs, trends['Transplant Rate'], 's-', color='green', label='Transplant Rate %')
    ax1.set_ylabel('Rate (%)', color='steelblue')
    ax1.tick_params(axis='y', labelcolor='steelblue')

    ax2 = ax1.twinx()
    ax2.plot(yrs, trends['Median OS'], 'D-', color='darkred', label='Median OS (months)')
    ax2.set_ylabel('Median OS (months)', color='darkred')
    ax2.tick_params(axis='y', labelcolor='darkred')

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1+lines2, labels1+labels2, frameon=False)
    ax1.set_title('Temporal Trends — Elderly Hepatobiliary Cancer (2000-2022)', fontweight='bold')
    _sax(ax1)
    _sax(ax2)

    fig.savefig('03_Analysis/figures/FigS4_TemporalTrends.png', dpi=300, bbox_inches='tight')
    fig.set_size_inches(W, fig.get_size_inches()[1] * W / fig.get_size_inches()[0])
    fig.savefig(os.path.join(FIGDIR, 'FigS4_TemporalTrends.png'), dpi=DPI, bbox_inches='tight', facecolor='white')
    fig.savefig(os.path.join(FIGDIR, 'FigS4_TemporalTrends.pdf'), bbox_inches='tight', facecolor='white')
    Image.open(os.path.join(FIGDIR, 'FigS4_TemporalTrends.png')).convert('RGB').save(os.path.join(FIGDIR, 'FigS4_TemporalTrends.tiff'), 'TIFF', compression='tiff_lzw', dpi=(DPI,DPI))
    plt.close()
    p("✓ Temporal trends saved")

    # Insight D: The "Segmental vs Larger" paradox across cancer types
    p("\n### D. Segmental vs Larger Resection: HCC vs ICC\n")
    for ct in ['HCC','ICC']:
        sub = df[df['cancer_type']==ct]
        p(f"\n**{ct}**")
        for s in ['Segmental_Resection','Larger_Resection']:
            grp = sub[sub['surgery_type']==s]
            p(f"  {s}: N={len(grp)}, median OS={grp['surv_months'].median():.0f}m, "
              f"mortality={grp['vital_dead'].mean()*100:.1f}%, "
              f"mean age={grp['age'].mean():.1f}")

        # HR comparison
        from lifelines import CoxPHFitter
        for s in ['Segmental_Resection','Larger_Resection']:
            sub2 = sub[sub['surgery_type'].isin(['None',s])]
            try:
                cph = CoxPHFitter()
                cph.fit(sub2[['surv_months','css_dead','surgery_any','age_c','stage_2','stage_3','stage_4',
                               'grade_poor','male','chemotherapy']],
                        'surv_months','css_dead')
                hr = np.exp(cph.params_['surgery_any'])
                ci = cph.confidence_intervals_.loc['surgery_any']
                ci_l, ci_h = np.exp(ci.iloc[0]), np.exp(ci.iloc[1])
                p(f"  {s} adjusted CSS HR: {hr:.2f} [{ci_l:.2f}-{ci_h:.2f}]")
            except Exception as e: p(f"  {s}: HR failed - {e}")

    # Insight E: Interaction effects
    p("\n### E. Significant Interactions (Cox with interactions)\n")
    try:
        int_df = df[['surv_months','css_dead','surgery_any','age_c','stage_2','stage_3','stage_4',
                     'is_icc','grade_poor','male','chemotherapy']].dropna()
        int_df['surg_x_age'] = int_df['surgery_any'] * int_df['age_c']
        int_df['surg_x_icc'] = int_df['surgery_any'] * int_df['is_icc']
        int_df['surg_x_stage4'] = int_df['surgery_any'] * int_df['stage_4']

        cph = CoxPHFitter(penalizer=0.01)
        cols = ['surgery_any','age_c','stage_2','stage_3','stage_4',
                'is_icc','grade_poor','male','chemotherapy',
                'surg_x_age','surg_x_icc','surg_x_stage4']
        cph.fit(int_df[['surv_months','css_dead']+cols], 'surv_months','css_dead')

        for c in cols:
            if c in cph.params_.index:
                hr = np.exp(cph.params_[c])
                pv = cph.summary.loc[c,'p']
                ps = "<0.001" if pv < 0.001 else f"{pv:.3f}"
                p(f"  {c}: HR={hr:.2f}, p={ps}")

        p(f"\n  C-index with interactions: {cph.concordance_index_:.3f}")
    except Exception as e:
        p(f"Interaction analysis failed: {e}")

print("\n=== ADVANCED ANALYSES COMPLETE ===")
