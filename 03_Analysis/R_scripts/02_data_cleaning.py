"""Phase 2: Data Cleaning & Feature Engineering for Hepatobiliary Cancer SEER Analysis

Reference: Zhang et al. (2020) Front Oncol + gastric cancer project pipeline
Enhanced with broader cancer types (HCC + ICC) and innovative features
"""
import pandas as pd
import numpy as np
import os, re, warnings
warnings.filterwarnings('ignore')

os.chdir(r'D:\Researching\SEER\hepatobiliary cancer')

# ============================================================
# 1. LOAD RAW DATA
# ============================================================
df = pd.read_csv(r'02_Data\raw\hepatobiliary cancer.csv', low_memory=False)
initial_n = len(df)
print(f"Raw data: {initial_n} rows")

# ============================================================
# 2. FILTER: Hepatobiliary cancers only (liver + intrahepatic bile duct)
# ============================================================
df = df[df['Primary Site - labeled'].isin(['C22.0-Liver', 'C22.1-Intrahepatic bile duct'])].copy()
print(f"After site filter: {len(df)} rows (dropped {initial_n - len(df)})")

# ============================================================
# 3. PARSE AGE
# ============================================================
def parse_age(x):
    if pd.isna(x): return np.nan
    m = re.search(r'(\d+)', str(x))
    return int(m.group(1)) if m else np.nan

df['age'] = df['Age recode with single ages and 90+'].apply(parse_age)
print(f"Age parsed: range {df['age'].min()}-{df['age'].max()}, NaN={df['age'].isna().sum()}")

# ============================================================
# 4. EXCLUDE: age < 18 (pediatric), missing age
# ============================================================
df = df[df['age'] >= 18].copy()
print(f"After age>=18: {len(df)}")

# ============================================================
# 5. HISTOLOGY: Classify HCC vs ICC vs Other
# ============================================================
hcc_codes = [8170, 8171, 8172, 8173, 8174, 8175, 8180]
icc_codes = [8160, 8161, 8162]

df['hist_code'] = pd.to_numeric(df['Histologic Type ICD-O-3'], errors='coerce')
df['cancer_type'] = 'Other'
df.loc[df['hist_code'].isin(hcc_codes), 'cancer_type'] = 'HCC'
df.loc[df['hist_code'].isin(icc_codes), 'cancer_type'] = 'ICC'
print(f"\nCancer type distribution:")
print(df['cancer_type'].value_counts())

# ============================================================
# 6. SURVIVAL OUTCOMES
# ============================================================
df['surv_months'] = pd.to_numeric(df['Survival months'], errors='coerce')
df['vital_dead'] = (df['Vital status recode (study cutoff used)'] == 'Dead').astype(int)
df['css_dead'] = (df['SEER cause-specific death classification'] == 'Dead (attributable to this cancer dx)').astype(int)

print(f"\nOS events: {df['vital_dead'].sum()} ({df['vital_dead'].mean()*100:.1f}%)")
print(f"CSS events: {df['css_dead'].sum()} ({df['css_dead'].mean()*100:.1f}%)")
print(f"Median survival: {df['surv_months'].median():.1f} months")

# ============================================================
# 7. UNIFIED STAGING: Merge across AJCC eras
# ============================================================
def unify_stage(row):
    """Merge AJCC 3rd/6th/7th/EOD 2018 into unified I-II-III-IV"""
    # Priority: 6th ed (2004-2015) > 7th ed (2010-2015) > EOD (2018+) > 3rd ed (1988-2003)
    for col in ['Derived AJCC Stage Group, 6th ed (2004-2015)',
                'Derived AJCC Stage Group, 7th ed (2010-2015)',
                'Derived EOD 2018 Stage Group Recode (2018+)']:
        val = str(row[col]).strip()
        # AJCC 6th/7th: I, IA, IB, II, IIA, IIB, III, IIIA, IIIB, IIIC, IIINOS, IV, IVA, IVB
        if val.startswith('III'):
            return 3
        elif val.startswith('II'):
            return 2
        elif val.startswith('I'):
            return 1
        elif val in ['IV', 'IVA', 'IVB'] or val.startswith('4'):
            return 4
        # EOD 2018 numeric: 1/1A/1B, 2, 3/3A/3B, 4/4A/4B
        try:
            s = int(float(val[0]))
            if 1 <= s <= 4:
                return s
        except:
            pass

    # Fall back to summary stage
    summary = str(row['Combined Summary Stage with Expanded Regional Codes (2004+)'])
    if 'Localized' in summary:
        return 1
    elif 'Regional' in summary:
        return 2
    elif 'Distant' in summary:
        return 4

    # 3rd edition fallback
    ajcc3 = str(row['AJCC stage 3rd edition (1988-2003)'])
    try:
        s = int(float(ajcc3))
        return min(max(s, 1), 4)
    except:
        pass

    return np.nan

df['stage'] = df.apply(unify_stage, axis=1)
print(f"\nUnified stage distribution:")
print(df['stage'].value_counts().sort_index())

# ============================================================
# 8. SURGERY CLASSIFICATION
# ============================================================
def classify_surgery(code):
    """Map SEER surgery codes to categories (paper: None/LD/SR/LR/LT)"""
    code = int(code) if not pd.isna(code) else 0
    if code == 0:
        return 'None'
    elif 10 <= code <= 27:
        return 'Local_Destruction'  # RFA, PEI, cryoablation
    elif 30 <= code <= 38:
        return 'Segmental_Resection'  # wedge/segmental
    elif 50 <= code <= 59:
        return 'Larger_Resection'  # lobectomy, extended lobectomy
    elif 60 <= code <= 66:
        return 'Transplant'
    elif code == 90:
        return 'Surgery_NOS'
    else:
        return 'Other_Surgery'

df['surgery_type'] = df['RX Summ--Surg Prim Site (1998+)'].apply(classify_surgery)
df['surgery_any'] = (df['surgery_type'] != 'None').astype(int)

print(f"\nSurgery type distribution:")
print(df['surgery_type'].value_counts())
print(f"Surgery any: {df['surgery_any'].sum()} ({df['surgery_any'].mean()*100:.1f}%)")

# ============================================================
# 9. DEMOGRAPHIC & CLINICAL FEATURES
# ============================================================

# Sex
df['male'] = (df['Sex'] == 'Male').astype(int)

# Race
def recode_race(x):
    x = str(x)
    if 'White' in x: return 'NHW'
    if 'Black' in x: return 'NHB'
    if 'Asian' in x or 'Pacific' in x: return 'NHAPI'
    if 'American Indian' in x or 'Alaska' in x: return 'NHAIAN'
    if 'Hispanic' in x: return 'Hispanic'
    return 'Other'

df['race'] = df['Race and origin recode (NHW, NHB, NHAIAN, NHAPI, Hispanic)'].apply(recode_race)

# Marital status
def recode_marital(x):
    x = str(x)
    if 'Married' in x: return 'Married'
    if 'Single' in x: return 'Single'
    if 'Divorced' in x: return 'Divorced'
    if 'Separated' in x: return 'Separated'
    if 'Widowed' in x: return 'Widowed'
    return 'Other'

df['marital'] = df['Marital status at diagnosis'].apply(recode_marital)
df['married'] = (df['marital'] == 'Married').astype(int)

# Grade (merge across eras)
def unify_grade(row):
    g1 = str(row['Grade Recode (thru 2017)'])
    g2 = str(row['Derived Summary Grade 2018 (2018+)'])
    g3 = str(row['Grade Clinical (2018+)'])
    g4 = str(row['Grade Pathological (2018+)'])

    for g in [g1, g2, g3, g4]:
        g = g.strip()
        if 'I' in g and 'II' not in g and not g.startswith('B'):
            return 1
        if 'II' in g and 'III' not in g:
            return 2
        if 'III' in g:
            return 3
        if 'IV' in g or 'Undifferentiated' in g or 'anaplastic' in g:
            return 4
        try:
            val = int(float(g))
            if 1 <= val <= 9:
                return val
        except:
            pass
    return np.nan

df['grade'] = df.apply(unify_grade, axis=1)
df['grade_poor'] = (df['grade'] >= 3).astype(int)

# Tumor size
df['tumor_size'] = pd.to_numeric(df['CS tumor size (2004-2015)'], errors='coerce')
df['tumor_size'] = df['tumor_size'].replace(999, np.nan)  # 999 = unknown

# === Tumor Size Summary (2016+) ===
# Parse numeric from Tumor Size Summary
def parse_size_2018(x):
    if pd.isna(x) or str(x).strip() == 'Blank(s)':
        return np.nan
    m = re.search(r'(\d+)', str(x))
    return int(m.group(1)) if m else np.nan

ts_2018 = df['Tumor Size Summary (2016+)'].apply(parse_size_2018)
df['tumor_size'] = df['tumor_size'].fillna(ts_2018)

# Nodes
df['nodes_examined'] = pd.to_numeric(df['Regional nodes examined (1988+)'], errors='coerce')
df['nodes_examined'] = df['nodes_examined'].replace(98, np.nan)
df['nodes_positive'] = pd.to_numeric(df['Regional nodes positive (1988+)'], errors='coerce')
df['nodes_positive'] = df['nodes_positive'].replace(98, np.nan)

# AFP
def recode_afp(x):
    x = str(x)
    if 'Positive' in x: return 'Positive'
    if 'Negative' in x: return 'Negative'
    if 'Borderline' in x: return 'Borderline'
    return 'Unknown'

df['afp'] = df['AFP Pretreatment Interpretation Recode (2010+)'].apply(recode_afp)

# Fibrosis
def recode_fibrosis(x):
    x = str(x)
    if '0-4' in x: return 'None-Moderate'
    if '5-6' in x or 'Advanced' in x or 'Cirrhosis' in x: return 'Advanced-Cirrhosis'
    return 'Unknown'

df['fibrosis'] = df['Fibrosis Score Recode (2010+)'].apply(recode_fibrosis)
df['cirrhosis'] = (df['fibrosis'] == 'Advanced-Cirrhosis').astype(int)

# Chemo / Radiation
df['chemotherapy'] = (df['Chemotherapy recode (yes, no/unk)'] == 'Yes').astype(int)
df['radiation'] = (df['Radiation recode'].str.contains('Beam|Radioisotope|implants|Combination', na=False)).astype(int)

# Year / income
df['year'] = df['Year of diagnosis']

def parse_income(x):
    x = str(x)
    m = re.search(r'\$?(\d+[,]?\d*)', x)
    if m:
        val = m.group(1).replace(',', '')
        return float(val) / 10000  # in $10K
    return np.nan

df['income_10k'] = df['Median household income inflation adj to 2023'].apply(parse_income)

# Age groups
df['age_group'] = pd.cut(df['age'], bins=[0, 50, 60, 65, 70, 75, 80, 100],
                         labels=['<50', '50-59', '60-64', '65-69', '70-74', '75-79', '80+'])

df['elderly'] = (df['age'] >= 65).astype(int)
df['very_elderly'] = (df['age'] >= 75).astype(int)

# ============================================================
# 10. BUILD ML DATASET (elderly patients, key features)
# ============================================================
# Filter: age >= 65 (elderly, matching paper), known stage, known survival
ml_df = df[(df['age'] >= 65) &
           (df['stage'].notna()) &
           (df['surv_months'].notna()) &
           (df['surv_months'] >= 0)].copy()

print(f"\n=== ML Dataset (age>=65, known stage/survival) ===")
print(f"N = {len(ml_df)}")

# Engineer ML features (matching gastric cancer project pattern)
ml_df['age_c'] = ml_df['age'] - 67  # centered at 67 (matching gastric project)

# Race dummies
for r in ['NHB', 'NHAPI', 'Hispanic', 'NHAIAN']:
    ml_df[f'race_{r.lower()}'] = (ml_df['race'] == r).astype(int)

# Stage dummies
for s in [2, 3, 4]:
    ml_df[f'stage_{s}'] = (ml_df['stage'] == s).astype(int)

# Surgery dummies
for st in ['Local_Destruction', 'Segmental_Resection', 'Larger_Resection', 'Transplant']:
    ml_df[f'surg_{st.lower()}'] = (ml_df['surgery_type'] == st).astype(int)

# Cancer type
ml_df['is_icc'] = (ml_df['cancer_type'] == 'ICC').astype(int)

# Feature list (matching gastric pattern + hepatobiliary-specific)
feature_cols = [
    'age_c', 'male', 'married',
    'race_nhb', 'race_nhapi', 'race_hispanic',
    'stage_2', 'stage_3', 'stage_4',
    'grade_poor',
    'is_icc',
    'tumor_size',
    'year',
    'income_10k',
    'surgery_any',
    'surg_local_destruction', 'surg_segmental_resection',
    'surg_larger_resection', 'surg_transplant',
    'chemotherapy', 'radiation',
    'cirrhosis',
]

# Outcome columns
outcome_cols = ['surv_months', 'vital_dead', 'css_dead']

# Select available columns
available_features = [c for c in feature_cols if c in ml_df.columns]
all_cols = available_features + outcome_cols + ['age', 'surgery_type', 'cancer_type', 'stage', 'race', 'grade', 'afp']

# Drop rows with NA in survival
ml_df = ml_df.dropna(subset=['surv_months', 'vital_dead'])

# Save
save_cols = [c for c in all_cols if c in ml_df.columns]
ml_df[save_cols].to_csv(r'02_Data\cleaned\hepatobiliary_elderly_clean.csv', index=False)
print(f"Cleaned data saved: {len(ml_df)} elderly patients, {len(available_features)} features")

# Feature summary
print(f"\nFeature list ({len(available_features)}):")
for c in available_features:
    missing = ml_df[c].isna().mean()
    print(f"  {c}: missing={missing:.1%}, mean={ml_df[c].mean():.3f}" if ml_df[c].dtype in ['float64','int64','int32','bool'] else f"  {c}: missing={missing:.1%}")

# Save full dataset (all ages) for comparison
df.to_csv(r'02_Data\cleaned\hepatobiliary_all_clean.csv', index=False)
print(f"\nFull dataset saved: {len(df)} patients")

print("\n=== DATA CLEANING COMPLETE ===")
