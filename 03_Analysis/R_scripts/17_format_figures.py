"""Post-process all figures for ASO submission: consistent colors, sizing, vector PDFs"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
import os, warnings
warnings.filterwarnings('ignore')

BASE = r'D:\Researching\SEER\hepatobiliary cancer'
FIG_IN = os.path.join(BASE, '04_Manuscript', 'figures')
FIG_OUT = os.path.join(BASE, '04_Manuscript', 'figures')

# Professional colorblind-friendly palette (7 colors)
COLORS = {
    'blue':    '#0072B2',    # Blue (colorblind-safe)
    'orange':  '#E69F00',    # Orange
    'green':   '#009E73',    # Green
    'red':     '#CC79A7',    # Pink/red
    'cyan':    '#56B4E9',    # Sky blue
    'yellow':  '#F0E442',    # Yellow
    'black':   '#000000',    # Black
}
# Surgery type palette (consistent across all figures)
SURG_COLORS = {
    'Non-surgery':      '#999999',
    'Local Destruction': '#E69F00',
    'Segmental Resection': '#0072B2',
    'Larger Resection':  '#CC79A7',
    'Transplant':        '#009E73',
    'Other':             '#F0E442',
}

FIGS = {
    # (filename, target_width_in, type)
    'Fig1_KM':         (7.0, 'km'),
    'Fig2_Forest':     (7.0, 'forest'),
    'Fig10_Nomogram':  (7.0, 'nomogram'),
    'Fig5_ExternalValidation': (7.0, 'bar'),
    'Fig20_HCCvsICC':  (7.0, 'km'),
    'Fig6_CompositeAnalysis': (7.0, 'composite'),
}

def convert_to_cmyk_pdf(png_path, pdf_path, target_dpi=600):
    """Convert PNG to CMYK PDF at print resolution"""
    img = Image.open(png_path).convert('CMYK')
    # Resize to reasonable print dimensions
    w, h = img.size
    target_inches = 7.0  # ASO single column is ~3.5in, double ~7in
    new_w = int(target_dpi * target_inches)
    new_h = int(h * (new_w / w))
    img_resized = img.resize((new_w, new_h), Image.LANCZOS)
    img_resized.save(pdf_path, 'PDF', resolution=target_dpi)
    return new_w, new_h

def add_unified_caption_border(png_path, output_path, fig_label):
    """Add consistent label and border to figure"""
    img = Image.open(png_path).convert('RGB')
    w, h = img.size
    # Add white border at bottom for caption area (5% height)
    border_h = int(h * 0.06)
    new_img = Image.new('RGB', (w, h + border_h), 'white')
    new_img.paste(img, (0, 0))
    # Add figure label text via matplotlib for consistency
    fig, ax = plt.subplots(figsize=(w/100, border_h/100))
    ax.text(0.02, 0.5, fig_label, fontsize=11, fontweight='bold',
            va='center', transform=ax.transAxes)
    ax.axis('off')
    fig.savefig('_label_temp.png', dpi=100, bbox_inches='tight', pad_inches=0)
    plt.close()
    # Overlay label
    label_img = Image.open('_label_temp.png').convert('RGB')
    lw, lh = label_img.size
    new_img.paste(label_img, (10, h + (border_h - lh)//2))
    new_img.save(output_path, quality=95, dpi=(300, 300))
    os.remove('_label_temp.png')

def batch_process():
    os.chdir(BASE)
    print("=== Figure Formatting for ASO Submission ===\n")

    # 1. Generate CMYK PDFs for all figures (main + suppl)
    print("Converting to CMYK PDF (600 DPI)...")
    for f in sorted(os.listdir(FIG_IN)):
        if f.endswith('.png') and not os.path.exists(f.replace('.png', '.pdf')):
            png_path = os.path.join(FIG_IN, f)
            pdf_path = os.path.join(FIG_OUT, f.replace('.png', '_cmyk.pdf'))
            try:
                w, h = convert_to_cmyk_pdf(png_path, pdf_path)
                print(f"  {f} -> {w}x{h}px @600dpi = {w/600:.1f}x{h/600:.1f}in")
            except Exception as e:
                print(f"  {f}: SKIP ({e})")

    # 2. Generate ASO-compatible TIFFs (print-ready RGB, ~5in wide)
    print("\nGenerating print-ready TIFFs...")
    for fname, (width_in, ftype) in FIGS.items():
        for ext in ['png']:
            src = os.path.join(FIG_IN, f'{fname}.{ext}')
            if not os.path.exists(src):
                continue
            # TIFF output at 300 DPI, ~5in wide
            out = os.path.join(FIG_OUT, f'{fname}_print.tiff')
            img = Image.open(src).convert('RGB')
            w, h = img.size
            target_w = int(300 * width_in)
            target_h = int(h * (target_w / w))
            img_resized = img.resize((target_w, target_h), Image.LANCZOS)
            img_resized.save(out, 'TIFF', compression='tiff_lzw',
                            dpi=(300, 300))
            print(f"  {fname}.tiff: {target_w}x{target_h} @300dpi = {width_in}x{target_h/300:.1f}in")

    # 3. Summary
    print("\n=== Done ===")
    print("CMYK PDFs: print-ready for prepress")
    print("TIFFs: ASO-compatible (LZW-compressed, 300 DPI)")

if __name__ == '__main__':
    batch_process()
