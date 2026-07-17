"""JAMA Surgery color palette — centralized source of truth for all figures.
Based on ggsci pal_jama + JAMA Network brand guidelines.
All 7 colors are colorblind-safe and grayscale-distinguishable."""

# ggsci JAMA palette: dark teal, orange, blue, red, green, purple, brown
C = {
    'non_surg':  '#374E55',  # dark teal/gray — neutral
    'local':     '#DF8F44',  # JAMA orange — primary accent
    'seg':       '#00A1D5',  # JAMA blue — complementary
    'larger':    '#B24745',  # JAMA red — accent
    'txp':       '#79AF97',  # JAMA green — positive outcome
    'other':     '#80796B',  # JAMA taupe/brown
    'hcc':       '#00A1D5',  # JAMA blue
    'icc':       '#6A6599',  # JAMA purple
    'ref':       '#80796B',  # JAMA taupe
    'head':      '#313131',  # JAMA mine shaft (brand dark gray)
}

SURG = {
    'None':               C['non_surg'],
    'Non-surgery':        C['non_surg'],
    'Local Destruction':  C['local'],
    'Segmental Resection': C['seg'],
    'Larger Resection':   C['larger'],
    'Transplant':         C['txp'],
    'Other':              C['other'],
}

SURG_SHORT = {
    'None':        C['non_surg'],
    'Local':       C['local'],
    'Segmental':   C['seg'],
    'Larger':      C['larger'],
    'Transplant':  C['txp'],
}

STAGE_C = {
    'I':   '#DF8F44',
    'II':  '#00A1D5',
    'III': '#B24745',
    'IV':  '#374E55',
}

AGE_C = {
    '65-69': '#DF8F44',
    '70-74': '#B24745',
    '75-79': '#6A6599',
    '80+':   '#374E55',
}

# prop_cycle: JAMA colors in logical order
PROP_CYCLE = ['#00A1D5','#DF8F44','#79AF97','#B24745','#6A6599','#80796B','#374E55']

JAMA_RC = {
    'font.family': 'sans-serif',
    'font.size': 8,
    'axes.titlesize': 9,
    'axes.labelsize': 8,
    'xtick.labelsize': 7.5,
    'ytick.labelsize': 7.5,
    'legend.fontsize': 7,
    'figure.dpi': 300,
    'axes.linewidth': 0.5,
    'xtick.major.width': 0.4,
    'ytick.major.width': 0.4,
    'axes.edgecolor': '#313131',
    'text.color': '#313131',
    'axes.labelcolor': '#313131',
    'xtick.color': '#374E55',
    'ytick.color': '#374E55',
}

# Color mapping: old hex → JAMA hex for batch replacement
COLOR_MAP = {
    '#7f8c8d': C['non_surg'],
    '#e67e22': C['local'],
    '#2980b9': C['seg'],
    '#3498db': C['seg'],
    '#8e44ad': C['icc'],
    '#27ae60': C['txp'],
    '#2c3e50': C['non_surg'],
    '#c0392b': C['larger'],
    '#d35400': C['local'],
    '#922b21': C['larger'],
    '#1a5276': C['non_surg'],
    '#1f6dad': C['other'],
    '#2b7bba': C['seg'],
    '#1a1a1a': C['head'],
    '#0072B2': C['seg'],
    '#E69F00': C['local'],
    '#009E73': C['txp'],
    '#CC79A7': C['larger'],
    '#999999': C['other'],
    '#56B4E9': C['seg'],
    '#F0E442': C['local'],
    '#bdc3c7': C['ref'],
    '#bdc3c7': C['ref'],
    '#222':    C['head'],
    '#555':    C['non_surg'],
    '#888':    C['other'],
    '#999':    C['other'],
    '#333':    C['non_surg'],
    '#1a6b3d': C['txp'],
    '#d4e6f1': '#E8F4F8',
}
