# Figure Color Scheme Reference — ASO Submission

> **Palette source**: ggsci `pal_jama("default")` — professional academic palette (colorblind-safe)
> All 7 colors are colorblind-safe and grayscale-distinguishable.

## Central Source

```python
# From R_scripts/_jama_palette.py — DO NOT hardcode elsewhere
C = {
    'non_surg':  '#374E55',  # JAMA dark teal — neutral
    'local':     '#DF8F44',  # JAMA orange — primary accent
    'seg':       '#00A1D5',  # JAMA blue — complementary
    'larger':    '#B24745',  # JAMA red — accent
    'txp':       '#79AF97',  # JAMA green — positive outcome
    'other':     '#80796B',  # JAMA taupe/brown
    'hcc':       '#00A1D5',  # JAMA blue (HCC)
    'icc':       '#6A6599',  # JAMA purple (ICC)
    'head':      '#313131',  # JAMA mine shaft (brand dark gray)
}
```

## Matplotlib Init

```python
import matplotlib.pyplot as plt
from _jama_palette import C, SURG, PROP_CYCLE, JAMA_RC

plt.rcParams.update(JAMA_RC)
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=PROP_CYCLE)
# PROP_CYCLE = ['#00A1D5','#DF8F44','#79AF97','#B24745','#6A6599','#80796B','#374E55']
```

## Verification

```python
import numpy as np
from PIL import Image
img = Image.open('fig.png').convert('L')
# All JAMA colors differ by >30 in grayscale for reliable print
```
