"""Final formatting: TIFF RGB 300dpi + resize to 174mm for ASO"""
import os
from PIL import Image

BASE = r'D:\Researching\SEER\hepatobiliary cancer'
FIG_DIR = os.path.join(BASE, '04_Manuscript', 'figures')
TARGET_INCH = 6.85  # 174mm
TARGET_DPI = 300

files = sorted(os.listdir(FIG_DIR))
png_files = [f for f in files if f.endswith('.png') and not f.endswith('_backup.png') and not '_cmyk' in f and not '_print' in f]

for f in png_files:
    path = os.path.join(FIG_DIR, f)
    img = Image.open(path)

    # 1. RGBA → RGB
    if img.mode == 'RGBA':
        # Composite over white background
        bg = Image.new('RGB', img.size, (255, 255, 255))
        bg.paste(img, mask=img.split()[3])  # use alpha as mask
        img = bg
    elif img.mode != 'RGB':
        img = img.convert('RGB')

    # 2. Resize to 174mm (6.85in) at 300dpi
    w, h = img.size
    new_w = int(TARGET_DPI * TARGET_INCH)
    new_h = int(h * (new_w / w))
    img_resized = img.resize((new_w, new_h), Image.LANCZOS)

    # 3. Save TIFF
    tiff_name = f.replace('.png', '.tiff')
    tiff_path = os.path.join(FIG_DIR, tiff_name)
    img_resized.save(tiff_path, 'TIFF', compression='tiff_lzw', dpi=(TARGET_DPI, TARGET_DPI))

    print(f"  {f}: {w}x{h} → {new_w}x{new_h}px @{TARGET_DPI}dpi = {TARGET_INCH:.2f}in")

print(f"\nDone: {len(png_files)} TIFF files at {TARGET_INCH}in width")
