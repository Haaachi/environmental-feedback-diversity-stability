"""
preview_palette.py
仅基于 annotation 文件预览真实数据的色板。
不需要测序/OD 数据。

用法
----
1. 把 Unified_species_annotations.xlsx (或 mapping .tsv/.csv) 放在
   ./data/ 目录下
2. 运行: python3 preview_palette.py
3. 输出: palette_preview.html (浏览器打开) + palette_preview.csv
"""

import os
import sys

# 确保能 import 主模块
BASE_DIR = os.environ.get("COMMUNITY_WORKSPACE", os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, BASE_DIR)

import pandas as pd
from process_abundance import (
    load_annotation_map,
    generate_color_palette,
)


def main():
    print("加载 annotation...")
    ann_df = load_annotation_map()

    m_taxa = set(ann_df[ann_df["experiment"] == "mortality"]["taxon"])
    t_taxa = set(ann_df[ann_df["experiment"] == "temperature"]["taxon"])

    print(f"  mortality:   {len(m_taxa)} taxa")
    print(f"  temperature: {len(t_taxa)} taxa")
    print(f"  总计:        {len(m_taxa) + len(t_taxa)} (两实验独立计数)")

    print("\n生成色板...")
    color_df = generate_color_palette(m_taxa, t_taxa, ann_df)

    # 写 CSV
    csv_path = os.path.join(BASE_DIR, "palette_preview.csv")
    color_df.to_csv(csv_path, index=False)
    print(f"  -> {csv_path}")

    # 写 HTML 预览
    html_path = os.path.join(BASE_DIR, "palette_preview.html")
    _write_html(color_df, html_path)
    print(f"  -> {html_path}")

    # 重复颜色检查
    for exp in ["mortality", "temperature"]:
        sub = color_df[color_df["experiment"] == exp]
        n_total = len(sub)
        n_unique = sub["color"].nunique()
        status = "OK" if n_total == n_unique else f"COLLISION ({n_total - n_unique} dup)"
        print(f"  [{exp}] {n_total} species, {n_unique} unique colors — {status}")


def _luma(hex_color):
    h = hex_color.lstrip("#")
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    return (r * 299 + g * 587 + b * 114) / 1000


def _write_html(color_df, out_path):
    rows_html = []

    for exp in ["mortality", "temperature"]:
        sub = color_df[color_df["experiment"] == exp].copy()
        if sub.empty:
            continue

        rows_html.append(
            f'<h2 style="font-size:18px; font-weight:500; margin:2rem 0 0.5rem;">'
            f'{exp} ({len(sub)} species)</h2>'
        )

        for phylum, sub_p in sub.groupby("phylum", sort=False):
            rows_html.append(
                f'<div style="font-size:13px; font-weight:500; color:#666; '
                f'margin:1rem 0 0.3rem; letter-spacing:0.02em;">{phylum}</div>'
            )

            for (family, genus), sub_g in sub_p.groupby(["family", "genus"], sort=False):
                cells = []
                for _, r in sub_g.iterrows():
                    text_col = "#000a" if _luma(r["color"]) > 140 else "#fffe"
                    cells.append(
                        f'<div title="{r["taxon"]} | {r["unified_annotation"]}" '
                        f'style="width:44px; height:24px; background:{r["color"]}; '
                        f'border-radius:3px; display:flex; align-items:center; '
                        f'justify-content:center; font-size:10px; color:{text_col}; '
                        f'font-weight:500;">{r["taxon"]}</div>'
                    )
                cells_html = "".join(cells)
                rows_html.append(
                    f'<div style="display:flex; align-items:center; gap:10px; '
                    f'margin-bottom:4px; font-size:13px;">'
                    f'<div style="width:320px; flex-shrink:0; color:#222;">'
                    f'<span style="color:#999;">{family}</span> · '
                    f'<span style="font-weight:500;">{genus}</span></div>'
                    f'<div style="display:flex; gap:3px; flex-wrap:wrap;">{cells_html}</div>'
                    f'</div>'
                )

    html = f"""<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8">
  <title>Palette preview</title>
  <style>
    body {{
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      max-width: 1200px;
      margin: 2rem auto;
      padding: 0 1rem;
      color: #222;
    }}
    h1 {{ font-size: 22px; font-weight: 500; }}
  </style>
</head>
<body>
  <h1>Palette preview</h1>
  <p style="color:#666; font-size:13px;">
    Hover a swatch to see its taxon ID. Each row is one genus.
  </p>
  {''.join(rows_html)}
</body>
</html>
"""

    with open(out_path, "w") as f:
        f.write(html)


if __name__ == "__main__":
    main()
