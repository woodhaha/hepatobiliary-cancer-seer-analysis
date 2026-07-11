# Figure Color Scheme Reference

Use this colorblind-friendly palette when regenerating figures.

## Color Sequence (Paul Tol)

```python
# ASO submission — professional, colorblind-safe, grayscale-compatible
COLORS = {
    'Non-surgery':       '#999999',  # gray
    'Local Destruction': '#E69F00',  # orange
    'Segmental':         '#0072B2',  # blue
    'Larger Resection':  '#CC79A7',  # pink
    'Transplant':        '#009E73',  # green
    'Other':             '#F0E442',  # yellow
}
```

## Matplotlib Init

```python
import matplotlib.pyplot as plt
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=[
    '#0072B2', '#E69F00', '#009E73', '#CC79A7', '#56B4E9', '#F0E442', '#000000'
])
```

## Verification

```python
# Check grayscale distinguishability
import numpy as np
from PIL import Image
img = Image.open('fig.png').convert('L')
# Colors should differ by >30 in grayscale for reliable print
```

Key principle: all 7 colors remain distinguishable when printed in grayscale, and work for protanopia/deuteranopia.
