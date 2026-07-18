"""Fig19_CONSORT — Clean rectangular boxes + thin arrows. ASO style."""
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os

BASE = r'D:\Researching\SEER\hepatobiliary cancer'
FIG = os.path.join(BASE, '04_Manuscript', 'figures')
os.makedirs(FIG, exist_ok=True)
W, DPI = 6.85, 300

plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 8, 'figure.dpi': DPI})

fig, ax = plt.subplots(figsize=(W, 5.5))
fig.patch.set_facecolor('white')
ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis('off')

def bx(x, y, w, h, txt, fs=7.5, fw='normal', tc='#222', fc='white', ec='#ccc'):
    ax.add_patch(mpatches.FancyBboxPatch((x-w/2, y-h/2), w, h,
        boxstyle="square", facecolor=fc, edgecolor=ec, lw=0.5, zorder=2))
    ax.text(x, y, txt, ha='center', va='center', fontsize=fs, fontweight=fw, color=tc, zorder=3)

def ar(y1, y2, x=0.5, c='#bbb'):
    ax.annotate('', xy=(x, y2), xytext=(x, y1),
                arrowprops=dict(arrowstyle='->', lw=0.4, color=c), zorder=1)

# Title
ax.text(0.5, 0.96, 'CONSORT Flow Diagram', fontsize=10, fontweight='bold', ha='center', color='#222')

# Main spine
bx(0.5, 0.86, 0.50, 0.060, 'SEER 18-Registry (2004–2022)\nN = 171,286', 7, 'normal', '#555')
ar(0.83, 0.73)

bx(0.5, 0.70, 0.50, 0.060, 'Patients with Complete Data\nN = 164,383', 7, 'normal', '#555')
ar(0.67, 0.57)

bx(0.5, 0.54, 0.50, 0.065, 'Elderly Study Cohort (Age ≥ 65)\nN = 76,110', 8.5, 'bold', '#1a1a1a', '#d4e6f1', '#2980b9')

# Exclusions — tight labels above arrows
ax.text(0.80, 0.785, '× 6,903 excluded', fontsize=5.5, color='#c0392b', ha='left')
ax.text(0.80, 0.635, '× 88,273 age < 65', fontsize=5.5, color='#c0392b', ha='left')

# Branch
ax.plot([0.27, 0.73], [0.50, 0.50], lw=0.3, color='#bbb')
for x in [0.27, 0.73]:
    ar(0.505, 0.42, x)
ax.text(0.27, 0.515, 'No surgery', fontsize=6, color='#555', ha='center')
ax.text(0.73, 0.515, 'Surgery', fontsize=6, color='#555', ha='center')

bx(0.27, 0.32, 0.34, 0.075, 'Non-Surgery\nn = 59,821 (78.6%)', 7.5, 'bold', '#c0392b')
bx(0.73, 0.32, 0.34, 0.075, 'Any Surgery\nn = 16,289 (21.4%)', 7.5, 'bold', '#1a6b3d')

ar(0.28, 0.18, 0.27)
ar(0.28, 0.18, 0.73)

bx(0.27, 0.10, 0.34, 0.065, 'Chemotherapy / Radiation\nBest Supportive Care', 6.5, 'normal', '#888')
bx(0.73, 0.08, 0.34, 0.095, 'Local Destruction: 10,771\nSegmental: 2,486  |  Larger: 699\nTransplant: 1,778', 6.5, 'normal', '#555')

name = 'Fig19_CONSORT'
plt.rcParams['tiff.compression'] = 'tiff_lzw'
fig.savefig(os.path.join(FIG, name+'.tiff'), dpi=DPI, bbox_inches='tight', facecolor='white')
fig.savefig(os.path.join(FIG, name+'.pdf'), bbox_inches='tight', facecolor='white')
plt.close()
print(f'  ✓ {name}')
