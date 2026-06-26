"""Phase 1: SEER Hepatobiliary Cancer Data Exploration"""
import pandas as pd
import numpy as np
import os, sys, warnings
warnings.filterwarnings('ignore')

os.chdir(r'D:\Researching\SEER\hepatobiliary cancer')
df = pd.read_csv(r'02_Data\raw\hepatobiliary cancer.csv', low_memory=False)

# Redirect output to file
out_path = r'03_Analysis\outputs\01_exploration_report.txt'
os.makedirs('03_Analysis/outputs', exist_ok=True)

with open(out_path, 'w', encoding='utf-8') as f:
    def p(*args, **kwargs):
        print(*args, **kwargs)
        print(*args, **kwargs, file=f)

    p(f"Shape: {df.shape}  |  Columns: {len(df.columns)}")
    p()

    # === Site Distribution ===
    p("=" * 60)
    p("SITE DISTRIBUTION (Primary Site - labeled)")
    p("=" * 60)
    site_counts = df['Primary Site - labeled'].value_counts()
    for site, cnt in site_counts.head(20).items():
        p(f"  {site}: {cnt} ({100*cnt/len(df):.1f}%)")
    p()

    # === Histology ===
    p("=" * 60)
    p("HISTOLOGIC TYPE")
    p("=" * 60)
    hist_counts = df['Histologic Type ICD-O-3'].value_counts()
    for h, cnt in hist_counts.head(15).items():
        p(f"  {h}: {cnt} ({100*cnt/len(df):.1f}%)")
    p()

    # === Year ===
    p("=" * 60)
    p("YEAR OF DIAGNOSIS")
    p("=" * 60)
    p(df['Year of diagnosis'].describe())
    p()
    yr = df['Year of diagnosis'].value_counts().sort_index()
    p(f"  Range: {yr.index.min()} - {yr.index.max()}")
    for y, cnt in yr.items():
        p(f"  {y}: {cnt}")
    p()

    # === Age ===
    p("=" * 60)
    p("AGE DISTRIBUTION")
    p("=" * 60)
    age = pd.to_numeric(df['Age recode with single ages and 90+'], errors='coerce')
    p(age.describe())
    p()
    p(f"  Age >= 65: {(age >= 65).sum()} ({(age >= 65).mean()*100:.1f}%)")
    p(f"  Age >= 70: {(age >= 70).sum()} ({(age >= 70).mean()*100:.1f}%)")
    p()

    # === Sex ===
    p("=" * 60)
    p("SEX")
    p("=" * 60)
    for v, c in df['Sex'].value_counts().items():
        p(f"  {v}: {c} ({100*c/len(df):.1f}%)")
    p()

    # === Race ===
    p("=" * 60)
    p("RACE/ETHNICITY")
    p("=" * 60)
    for v, c in df['Race and origin recode (NHW, NHB, NHAIAN, NHAPI, Hispanic)'].value_counts().items():
        p(f"  {v}: {c} ({100*c/len(df):.1f}%)")
    p()

    # === Surgery ===
    p("=" * 60)
    p("SURGERY (RX Summ--Surg Prim Site)")
    p("=" * 60)
    surg = df['RX Summ--Surg Prim Site (1998+)']
    p(surg.describe())
    p()
    surg_counts = surg.value_counts().sort_index()
    for v, c in surg_counts.head(30).items():
        p(f"  Code {v}: {c} ({100*c/len(df):.1f}%)")
    p()

    # === Stage ===
    p("=" * 60)
    p("AJCC STAGE 6th ed (2004-2015)")
    p("=" * 60)
    for v, c in df['Derived AJCC Stage Group, 6th ed (2004-2015)'].value_counts().items():
        p(f"  {v}: {c} ({100*c/len(df):.1f}%)")
    p()

    p("=" * 60)
    p("COMBINED SUMMARY STAGE")
    p("=" * 60)
    for v, c in df['Combined Summary Stage with Expanded Regional Codes (2004+)'].value_counts().head(15).items():
        p(f"  {v}: {c} ({100*c/len(df):.1f}%)")
    p()

    # === Survival ===
    p("=" * 60)
    p("SURVIVAL MONTHS")
    p("=" * 60)
    surv = pd.to_numeric(df['Survival months'], errors='coerce')
    p(surv.describe())
    p()

    # === Vital Status ===
    p("=" * 60)
    p("VITAL STATUS")
    p("=" * 60)
    for v, c in df['Vital status recode (study cutoff used)'].value_counts().items():
        p(f"  {v}: {c} ({100*c/len(df):.1f}%)")
    p()

    # === CSS ===
    p("=" * 60)
    p("SEER CAUSE-SPECIFIC DEATH")
    p("=" * 60)
    for v, c in df['SEER cause-specific death classification'].value_counts().items():
        p(f"  {v}: {c} ({100*c/len(df):.1f}%)")
    p()

    # === Grade ===
    p("=" * 60)
    p("GRADE")
    p("=" * 60)
    for v, c in df['Grade Recode (thru 2017)'].value_counts().items():
        p(f"  {v}: {c} ({100*c/len(df):.1f}%)")
    p()

    # === Tumor Size ===
    p("=" * 60)
    p("TUMOR SIZE (CS tumor size 2004-2015)")
    p("=" * 60)
    ts = pd.to_numeric(df['CS tumor size (2004-2015)'], errors='coerce')
    p(ts.describe())
    p(f"  999 (Unknown): {(ts==999).sum()}")
    p()

    # === AFP ===
    p("=" * 60)
    p("AFP")
    p("=" * 60)
    for v, c in df['AFP Pretreatment Interpretation Recode (2010+)'].value_counts().items():
        p(f"  {v}: {c} ({100*c/len(df):.1f}%)")
    p()

    # === Fibrosis ===
    p("=" * 60)
    p("FIBROSIS SCORE")
    p("=" * 60)
    for v, c in df['Fibrosis Score Recode (2010+)'].value_counts().items():
        p(f"  {v}: {c} ({100*c/len(df):.1f}%)")
    p()

    # === Chemotherapy ===
    p("=" * 60)
    p("CHEMOTHERAPY")
    p("=" * 60)
    for v, c in df['Chemotherapy recode (yes, no/unk)'].value_counts().items():
        p(f"  {v}: {c} ({100*c/len(df):.1f}%)")
    p()

    # === Radiation ===
    p("=" * 60)
    p("RADIATION")
    p("=" * 60)
    for v, c in df['Radiation recode'].value_counts().head(10).items():
        p(f"  {v}: {c} ({100*c/len(df):.1f}%)")
    p()

    # === Marital status ===
    p("=" * 60)
    p("MARITAL STATUS")
    p("=" * 60)
    for v, c in df['Marital status at diagnosis'].value_counts().items():
        p(f"  {v}: {c} ({100*c/len(df):.1f}%)")
    p()

    # === Income ===
    p("=" * 60)
    p("MEDIAN HOUSEHOLD INCOME")
    p("=" * 60)
    for v, c in df['Median household income inflation adj to 2023'].value_counts().items():
        p(f"  {v}: {c} ({100*c/len(df):.1f}%)")

print(f"DONE: {out_path}")
