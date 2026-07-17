# JAMA Surgery Submission Package

**Manuscript**: Comparison of Segmental Resection vs Major Hepatectomy in Elderly Hepatobiliary Cancer
**Date**: July 17, 2026

---

## Required Files

| Item | File | Notes |
|------|------|-------|
| **Cover Letter** | `cover_letter_jamasurg.md` | Addressed to Gerard Doherty, MD, Editor-in-Chief |
| **Manuscript** | `manuscript_jamasurg.docx` | Word format, 12pt Times New Roman, double-spaced, 2,847 words |
| **STROBE Checklist** | `STROBE_checklist.md` | 22/22 items complete (updated for JAMA version) |

### Main Figures (separate files, 1 per figure)

| Figure | File (preferred format) | Backup | Description |
|--------|------------------------|--------|-------------|
| Fig 1 | `Fig1_KM.tiff` (300 DPI, RGB, LZW) | `.png` / `.pdf` | Kaplan-Meier: CSS by surgery type |
| Fig 2 | `Fig2_Forest.tiff` (300 DPI, RGB, LZW) | `.png` / `.pdf` | Multivariable Cox forest plot |
| Fig 3 | `Fig5_HBI.tiff` (300 DPI, RGB, LZW) | `.png` / `.pdf` | HBI risk score + decision support |

All figures: 300 DPI, ~5 in wide, TIFF with LZW compression, RGB color.
All ≤500 KB per § file size limit.

### Online-Only Supplementary Material

**eFigures** (all available as .tiff + .png + .pdf):

| eFig | File | Description |
|------|------|-------------|
| eFig1 | `eFig1_Consort.tiff` | CONSORT patient selection diagram |
| eFig2 | `eFig2_SHAP.tiff` | SHAP feature importance analysis |
| eFig3 | `eFig3_ModelAnalysis.tiff` | 4-model comparison (Cox, RSF, XGBoost, DeepSurv) |
| eFig4 | `eFig4_Comprehensive.tiff` | Temporal trends + age-dependent surgery benefit |
| eFig5 | `eFig5_LandmarkFrailty.tiff` | Landmark conditional survival + frailty analysis |
| eFig6 | `eFig6_Calibration.tiff` | Model calibration at 12/36/60 months |
| eFig7 | `eFig7_Subgroup.tiff` | Expanded subgroup interaction forest plot |
| eFig8 | `eFig8_Overview.tiff` | Comprehensive study overview (8 panels) |

**eTables**: Refer to supplementary PDF (`supplementary_jama.pdf`) for full supplementary tables (eTable 1–10).

### Supplementary PDF
- `supplementary_jama.pdf` — Contains all supplementary evidence: full baseline table, complete Cox results, E-value stratification, model specification tests, leave-one-out analysis, geographic variation, Schoenfeld residuals, hyperparameters, decision matrix, included vs excluded comparison.

---

## JAMA Surgery Format Compliance Checklist

| Requirement | Status |
|-------------|--------|
| Title ≤100 characters (86 ✓) | ✅ |
| Key Points before Abstract (97 words ✓) | ✅ |
| Structured abstract ≤350 words (350 ✓) | ✅ |
| Body text ≤3000 words (1,926 ✓) | ✅ |
| ≤5 tables/figures total (2 tables + 3 figures = 5 ✓) | ✅ |
| References ≤50 (29 ✓) | ✅ |
| Data sharing statement included | ✅ |
| STROBE checklist complete (22/22) | ✅ |
| Figures as separate files (not embedded) | ✅ |
| Figures in TIFF ≥300 DPI, RGB | ✅ |
| Manuscript in Word .docx format | ✅ |
| Double-spaced, 1in margins, Times New Roman | ✅ |
| Cover letter addressed to JAMA Surgery editor | ✅ |
| No study type in title | ✅ |
| No declarative title | ✅ |
| TRIPOD Type 3 reporting | ✅ |

---

## Submission Strategy

**Primary target**: JAMA Surgery (Original Investigation)

**Fallback if rejected**:
1. JAMA Network Open (higher acceptance rate, receptive to SEER studies)
2. Annals of Surgery (surgical audience, higher bar for novelty)
3. HPB or JCO Oncology (specialty journal, likely acceptance)

**Next steps**:
1. Submit via JAMA Surgery online submission system (https://ms.manuscriptcentral.com/jamasurgery)
2. Await editorial decision (typical: 1–2 weeks for initial screening, 4–8 weeks for peer review)
3. If desk rejected, transfer materials to JAMA Network Open (same format)
