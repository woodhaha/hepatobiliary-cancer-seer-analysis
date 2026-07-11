"""Inject Paul Tol colorblind palette into all plotting scripts"""
import os, re

BASE = r'D:\Researching\SEER\hepatobiliary cancer\03_Analysis\R_scripts'
RC_LINE = "plt.rcParams['axes.prop_cycle'] = plt.cycler(color=['#0072B2','#E69F00','#009E73','#CC79A7','#56B4E9','#F0E442','#000000'])"

for f in sorted(os.listdir(BASE)):
    if not f.endswith('.py') or f.startswith('_') or f == '17_format_figures.py':
        continue
    fpath = os.path.join(BASE, f)
    with open(fpath, 'r', encoding='utf-8') as fh:
        content = fh.read()

    # Skip if already injected
    if 'axes.prop_cycle' in content:
        print(f"  {f}: already has color scheme, skip")
        continue

    # Find import matplotlib line and add rcParams after it
    lines = content.split('\n')
    new_lines = []
    injected = False
    for i, line in enumerate(lines):
        new_lines.append(line)
        if not injected and 'import matplotlib.pyplot as plt' in line and 'plt.rcParams' not in content:
            # Check next line isn't already rcParams
            if i + 1 >= len(lines) or 'axes.prop_cycle' not in lines[i + 1]:
                new_lines.append(RC_LINE)
                injected = True

    if injected:
        with open(fpath, 'w', encoding='utf-8') as fh:
            fh.write('\n'.join(new_lines))
        print(f"  {f}: injected color scheme")
    else:
        print(f"  {f}: no matplotlib import found, skip")

print("\nDone. Injected color scheme into figure scripts.")
