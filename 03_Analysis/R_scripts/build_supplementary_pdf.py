"""Build JAMA Network Open supplementary materials PDF with embedded figures."""
import os, sys
from PIL import Image
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import textwrap

BASE = r'D:\Researching\SEER\hepatobiliary cancer'
FIG = os.path.join(BASE, '04_Manuscript', 'figures')
OUT = os.path.join(BASE, '05_Submission')

# Maps eFigure number → PNG path → title
EFIGS = [
    (1,  'eFig1_Consort.png',           'eFigure 1. CONSORT Flow Diagram'),
    (2,  'eFig2_SHAP.png',              'eFigure 2. SHAP Feature Importance'),
    (3,  'eFig3_ModelAnalysis.png',      'eFigure 3. Model Analysis: LASSO Selection + DCA + Competing Risk'),
    (4,  'eFig4_Comprehensive.png',      'eFigure 4. Age-Dependent Surgery Benefit & Temporal Trends'),
    (5,  'eFig5_LandmarkFrailty.png',    'eFigure 5. Landmark Conditional Survival & Frailty Stratification'),
    (6,  'eFig6_Calibration.png',        'eFigure 6. Model Calibration: Brier Scores & Risk Stratification'),
    (7,  'eFig7_Subgroup.png',           'eFigure 7. Subgroup Interaction Forest Plot'),
    (8,  'eFig8_Overview.png',           'eFigure 8. Comprehensive Analysis Overview (8 panels)'),
]

W, H = 8.27, 11.69  # A4

def make_page(num, title, img_path=None, img_caption=None):
    """Create a matplotlib figure page for the supplementary PDF."""
    fig, ax = plt.subplots(figsize=(W, H))
    ax.set_xlim(0, W); ax.set_ylim(0, H)
    ax.axis('off')

    # Header
    ax.text(W/2, H-0.3, title, fontsize=12, fontweight='bold', ha='center', va='top')

    y = H - 0.8

    if img_path and os.path.exists(img_path):
        try:
            img = Image.open(img_path)
            # Scale to fit page width
            iw, ih = img.size
            scale = min((W-2) / iw, (H - 2.5) / ih)
            nw, nh = iw * scale, ih * scale
            ax.imshow(img, extent=[1, 1+nw, y-nh, y], aspect='auto')
            y -= nh + 0.3
        except Exception as e:
            ax.text(1, y-0.5, f'[Image could not be loaded: {e}]', fontsize=9, color='red')
            y -= 1
    else:
        ax.text(1, y-0.3, '[Refer to Supplementary Materials text for details]', fontsize=9, style='italic')

    if img_caption:
        y -= 0.1
        for line in textwrap.wrap(img_caption, width=100):
            ax.text(1, y, line, fontsize=8, color='#333')
            y -= 0.2

    return fig

print("Building supplementary PDF...")

# ─── Title / header page ───
figs = []
fig, ax = plt.subplots(figsize=(W, H))
ax.set_xlim(0, W); ax.set_ylim(0, H)
ax.axis('off')
ax.text(W/2, H-1, 'Supplementary Materials', fontsize=18, fontweight='bold', ha='center')
ax.text(W/2, H-1.5, 'Surgical Resection and Extent in Elderly Patients With Hepatobiliary Cancer',
        fontsize=11, ha='center')
ax.text(W/2, H-1.9, 'JAMA Network Open', fontsize=10, ha='center', color='#555')
ax.text(1, H-3, 'Contents:', fontsize=11, fontweight='bold')
contents = [
    'eFigures 1-8',
    '  eFigure 1: CONSORT Flow Diagram',
    '  eFigure 2: SHAP Feature Importance',
    '  eFigure 3: Model Analysis (LASSO + DCA + Competing Risk)',
    '  eFigure 4: Age-Dependent Surgery Benefit & Temporal Trends',
    '  eFigure 5: Landmark Analysis & Frailty Stratification',
    '  eFigure 6: Model Calibration',
    '  eFigure 7: Subgroup Interaction Forest',
    '  eFigure 8: Comprehensive Analysis Overview',
    '',
    'eTables 1-9 (included in supplementary_materials.md)',
    '  eTable 1: Full Baseline Characteristics',
    '  eTable 2: Complete Multivariate Cox Regression',
    '  eTable 3: Stratified E-values',
    '  eTable 4: Model Specification Robustness',
    '  eTable 5: Leave-One-Out Sensitivity Analysis',
    '  eTable 6: Geographic Practice Variation',
    '  eTable 7: Model Hyperparameters',
    '  eTable 8: Schoenfeld Residual Tests',
    '  eTable 9: Clinical Decision Matrix',
]
for i, line in enumerate(contents):
    ax.text(1.5, H-3.3 - i*0.25, line, fontsize=9, va='top')
figs.append(fig)

# ─── eFigure pages ───
for num, img_name, title in EFIGS:
    img_path = os.path.join(FIG, img_name) if img_name else None
    f = make_page(num, title, img_path=img_path, img_caption=None)
    figs.append(f)

# ─── STROBE Checklist page ───
fig, ax = plt.subplots(figsize=(W, H))
ax.set_xlim(0, W); ax.set_ylim(0, H)
ax.axis('off')
ax.text(W/2, H-0.3, 'STROBE Checklist — Cohort Study', fontsize=12, fontweight='bold', ha='center')
strobe = [
    ('Title & Abstract', [
        '1a: Indicate study design in title ✓',
        '1b: Informative abstract ✓',
    ]),
    ('Introduction', [
        '2: Scientific background ✓',
        '3: Specific objectives ✓',
    ]),
    ('Methods', [
        '4: Key study design elements ✓',
        '5: Setting, locations, dates ✓',
        '6a: Eligibility criteria ✓',
        '6b: Follow-up methods ✓',
        '7: Outcomes/exposures defined ✓',
        '8: Data sources ✓',
        '9: Bias mitigation (PSM, E-value) ✓',
        '10: Study size (CONSORT) ✓',
        '11: Quantitative variables ✓',
        '12a: Statistical methods ✓',
        '12b: Subgroup methods ✓',
        '12c: Missing data addressed ✓',
        '12d: Sensitivity analyses ✓',
    ]),
    ('Results', [
        '13a: Individuals at each stage (CONSORT) ✓',
        '13b: Reasons for exclusion ✓',
        '13c: Flow diagram ✓',
        '14a: Participant characteristics ✓',
        '14b: Missing data counts ✓',
        '14c: Follow-up time ✓',
        '15: Outcome events ✓',
        '16a: Unadjusted estimates ✓',
        '16b: Category boundaries ✓',
        '17: Other analyses (ML, external validation) ✓',
    ]),
    ('Discussion', [
        '18: Key results with objectives ✓',
        '19: Limitations ✓',
        '20: Cautious interpretation ✓',
        '21: Generalisability ✓',
    ]),
    ('Other', [
        '22: Funding source ✓',
    ]),
]
y = H - 0.8
for section, items in strobe:
    ax.text(1.5, y, section, fontsize=10, fontweight='bold'); y -= 0.3
    for item in items:
        ax.text(2, y, f'• {item}', fontsize=7.5); y -= 0.18
    y -= 0.05
ax.text(W/2, 0.5, 'Overall: 22/22 items addressed ✓', fontsize=10, fontweight='bold', ha='center')
figs.append(fig)

# ─── Supplementary References page ───
fig, ax = plt.subplots(figsize=(W, H))
ax.set_xlim(0, W); ax.set_ylim(0, H)
ax.axis('off')
ax.text(W/2, H-0.3, 'Supplementary References', fontsize=12, fontweight='bold', ha='center')
refs = [
    'S1. TCGA Research Network. Comprehensive and integrative genomic characterization of hepatocellular carcinoma. Cell. 2017;169(7):1327-1341.',
    'S2. Totoki Y, et al. Trans-ancestry mutational landscape of hepatocellular carcinoma genomes. Nat Genet. 2014;46:1267-1273.',
    'S3. VanderWeele TJ, Ding P. Sensitivity analysis in observational research: introducing the E-value. Ann Intern Med. 2017;167(4):268-274.',
    'S4. Baiocchi M, Cheng J, Small DS. Instrumental variable methods for causal inference. Stat Med. 2014;33(13):2297-2340.',
]
y = H - 0.8
for ref in refs:
    wrapped = textwrap.wrap(ref, width=110)
    for line in wrapped:
        ax.text(1.5, y, line, fontsize=8); y -= 0.18
    y -= 0.1
figs.append(fig)

# ─── Save as multipage PDF via PdfPages ───
outpath = os.path.join(OUT, 'supplementary_jama.pdf')
with PdfPages(outpath) as pdf:
    for f in figs:
        pdf.savefig(f, dpi=200, bbox_inches='tight', pad_inches=0.3, facecolor='white')
plt.close('all')
print(f'✅ Supplementary PDF saved: {outpath}')
print(f'   Pages: {len(figs)}')
