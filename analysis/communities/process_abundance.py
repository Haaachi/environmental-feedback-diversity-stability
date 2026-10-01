"""
process_abundance.py
主目录: /home/hachi/Tem_mortality_workspace

功能
----
1. 从统一测序 Excel 和 OD Excel 中解析数据
2. 计算 relative abundance 和 absolute abundance
3. 加载 species annotation
4. 输出 processed/{experiment}/W*/abs_abundance.csv
5. 输出 color_map.csv

颜色映射设计（保留最初版风格 + 轻量分类约束）
-----------------------------------------------
**风格**（完全沿用你最初版）：
  - HSL 色彩空间
  - 跳过红色 (0.00-0.06, 0.92-1.00) 和纯绿色 (0.28-0.42)
  - 6 种 sat/lum 组合循环，中等亮度
  - interleaving 排列：相邻 species 拿到 hue 距离最远的位置
  - mortality / temperature 用 hue_start offset 错开 (0.0 / 0.5)

**唯一的分类约束**：
  - species 在 hue 轴上的位置由"taxon 字典序"改成
    "phylum → family → genus → taxon 数字号"排序后的位次
  - 这样同 family/genus 的 species 拿到相邻 hue slot，色相相近
  - 但 sat/lum 组合循环 + interleaving 仍然让每个 species 清晰可分

**为什么不再做 family 加权 / phylum 分块 / L 梯度**：
  - 那些约束会把全局 hue 空间切碎，导致 species 多时颜色趋同
  - 你最初版的"全 species 平铺"之所以颜色丰富，正是因为没有这种切分
  - 这里只用 annotation 来决定"谁排在谁旁边"，不切空间
"""

import re
import os
import colorsys
import numpy as np
import pandas as pd
import openpyxl


BASE_DIR = os.environ.get("COMMUNITY_WORKSPACE", os.path.dirname(os.path.abspath(__file__)))
DATA_DIR  = os.path.join(BASE_DIR, "data")
OUT_DIR   = os.path.join(BASE_DIR, "processed")

ANNOTATION_CANDIDATES = [
    os.path.join(DATA_DIR, "Unified_species_annotations.xlsx"),
    os.path.join(DATA_DIR, "Unified_species_annotation_mapping.xlsx"),
    os.path.join(DATA_DIR, "Unified_species_annotation_mapping.tsv"),
    os.path.join(DATA_DIR, "Unified_species_annotation_mapping.csv"),
]

ANNOTATION_XLSX_SHEET = "Species_mapping"

SEQ_FILES = {
    "mortality":   os.path.join(DATA_DIR, "Mortality_Unified.xlsx"),
    "temperature": os.path.join(DATA_DIR, "Temperature_Unified.xlsx"),
}

OD_FILES = {
    "mortality":   os.path.join(DATA_DIR, "Mortality_OD_Unified.xlsx"),
    "temperature": os.path.join(DATA_DIR, "Temperature_OD_Unified.xlsx"),
}

OD_BACKGROUND = {
    "mortality":   0.040,
    "temperature": 0.035,
}

SEQ_COLS = {
    "D1(R1)":  (1, 1),  "D2(R1)":  (1, 2),  "D3(R1)":  (1, 3),
    "D4(R1)":  (1, 4),  "D5(R1)":  (1, 5),  "D6(R1)":  (1, 6),
    "D7(R1)":  (1, 7),  "D8(R1)":  (1, 8),  "D9(R1)":  (1, 9),
    "D10(R1)": (1, 10),
    "D8(R2)":  (2, 8),  "D9(R2)":  (2, 9),  "D10(R2)": (2, 10),
    "D8(R3)":  (3, 8),  "D9(R3)":  (3, 9),  "D10(R3)": (3, 10),
}


# =====================================================================
#  Sequencing & OD parsing
# =====================================================================

def parse_seq_excel(path, experiment):
    wb = openpyxl.load_workbook(path)
    records = []

    for sheet in wb.sheetnames:
        if sheet == "README":
            continue

        condition = sheet
        ws = wb[sheet]

        current_comm = None
        headers = None

        for row in ws.iter_rows(values_only=True):
            c0 = str(row[0]).strip() if row[0] is not None else ""

            if c0.startswith("Community"):
                m = re.search(r"Community\s+(\d+)", c0)
                current_comm = int(m.group(1)) if m else None
                headers = None
                continue

            if c0 in ("Taxon", "Species"):
                headers = list(row)
                continue

            if headers and current_comm and c0 and c0 not in ("Taxon", "Species"):
                for col_i, hdr in enumerate(headers):
                    if hdr in SEQ_COLS:
                        replica, day = SEQ_COLS[hdr]
                        val = row[col_i]

                        if val is None:
                            reads = np.nan
                        else:
                            try:
                                reads = int(float(val))
                            except Exception:
                                reads = 0

                        records.append({
                            "experiment": experiment,
                            "condition":  condition,
                            "community":  current_comm,
                            "taxon":      c0,
                            "replica":    replica,
                            "day":        day,
                            "reads":      reads,
                        })

    df = pd.DataFrame(records)

    if df.empty:
        return df

    expected = pd.DataFrame([
        {"replica": replica, "day": day}
        for replica, day in sorted(set(SEQ_COLS.values()))
    ])

    keys = (
        df[["experiment", "condition", "community", "taxon"]]
        .drop_duplicates()
    )

    full_index = keys.merge(expected, how="cross")

    df = full_index.merge(
        df,
        on=["experiment", "condition", "community", "taxon", "replica", "day"],
        how="left",
    )

    return df


def parse_od_excel(path, experiment):
    bg = OD_BACKGROUND[experiment]
    wb = openpyxl.load_workbook(path)
    records = []

    for sheet in wb.sheetnames:
        if sheet == "README":
            continue

        condition = sheet
        ws = wb[sheet]

        current_comm = None

        for row in ws.iter_rows(values_only=True):
            c0 = row[0]

            if (
                c0 is not None
                and str(c0).startswith("Community")
                and str(c0) != "Community/Replicate"
            ):
                m = re.search(r"Community(\d+)", str(c0))
                current_comm = int(m.group(1)) if m else None

                rep_label = str(row[12]).strip() if row[12] else ""

                if rep_label.startswith("R"):
                    m2 = re.search(r"\d+", rep_label)
                    if not m2:
                        continue

                    rep_num = int(m2.group())

                    for day_i, day in enumerate(range(1, 11), start=2):
                        od_val = row[day_i]
                        od = max(0.0, float(od_val) - bg) if od_val is not None else None

                        records.append({
                            "experiment": experiment,
                            "condition":  condition,
                            "community":  current_comm,
                            "replica":    rep_num,
                            "day":        day,
                            "OD":         round(od, 6) if od is not None else None,
                        })

                continue

            if c0 is None and current_comm is not None:
                rep_label = str(row[12]).strip() if row[12] else ""

                if rep_label.startswith("R"):
                    m2 = re.search(r"\d+", rep_label)
                    if not m2:
                        continue

                    rep_num = int(m2.group())

                    for day_i, day in enumerate(range(1, 11), start=2):
                        od_val = row[day_i]
                        od = max(0.0, float(od_val) - bg) if od_val is not None else None

                        records.append({
                            "experiment": experiment,
                            "condition":  condition,
                            "community":  current_comm,
                            "replica":    rep_num,
                            "day":        day,
                            "OD":         round(od, 6) if od is not None else None,
                        })

    return pd.DataFrame(records)


def compute_abundance(seq_df, od_df):
    totals = (
        seq_df.dropna(subset=["reads"])
        .groupby(["experiment", "condition", "community", "replica", "day"])["reads"]
        .sum()
        .reset_index()
        .rename(columns={"reads": "total_reads"})
    )

    df = seq_df.merge(
        totals,
        on=["experiment", "condition", "community", "replica", "day"],
        how="left",
    )

    df["rel_abund"] = np.where(
        df["reads"].isna() | (df["total_reads"] == 0),
        np.nan,
        df["reads"] / df["total_reads"] * 100,
    )

    df = df.merge(
        od_df,
        on=["experiment", "condition", "community", "replica", "day"],
        how="left",
    )

    df["abs_abund"] = np.where(
        df["rel_abund"].isna() | df["OD"].isna(),
        np.nan,
        df["rel_abund"] / 100.0 * df["OD"],
    )

    df["rel_abund"] = df["rel_abund"].round(6)
    df["abs_abund"] = df["abs_abund"].round(8)

    return df


# =====================================================================
#  Annotation handling
# =====================================================================

GENUS_TAXONOMY = {
    "Unclassified_Enterobacteriaceae": ("Proteobacteria", "Enterobacteriaceae"),
    "Klebsiella": ("Proteobacteria", "Enterobacteriaceae"),
    "Pantoea": ("Proteobacteria", "Erwiniaceae"),
    "Phytobacter": ("Proteobacteria", "Erwiniaceae"),

    "Bacillus": ("Firmicutes", "Bacillaceae"),
    "Lysinibacillus": ("Firmicutes", "Planococcaceae"),
    "Unclassified_Planococcaceae": ("Firmicutes", "Planococcaceae"),
    "Exiguobacterium": ("Firmicutes", "Exiguobacteraceae"),
    "Paenibacillus": ("Firmicutes", "Paenibacillaceae"),
    "Staphylococcus": ("Firmicutes", "Staphylococcaceae"),

    "Chryseobacterium": ("Bacteroidota", "Weeksellaceae"),
    "Sphingobacterium": ("Bacteroidota", "Sphingobacteriaceae"),
    "Flavobacterium": ("Bacteroidota", "Flavobacteriaceae"),

    "Curtobacterium": ("Actinobacteriota", "Microbacteriaceae"),
    "Microbacterium": ("Actinobacteriota", "Microbacteriaceae"),

    "Burkholderia": ("Proteobacteria", "Burkholderiaceae"),
    "Burkholderia-Caballeronia-Paraburkholderia": ("Proteobacteria", "Burkholderiaceae"),
    "Pandoraea": ("Proteobacteria", "Burkholderiaceae"),
    "Unclassified_Comamonadaceae": ("Proteobacteria", "Comamonadaceae"),
    "Aeromonas": ("Proteobacteria", "Aeromonadaceae"),
    "Pseudomonas": ("Proteobacteria", "Pseudomonadaceae"),
    "Acinetobacter": ("Proteobacteria", "Moraxellaceae"),
    "Stenotrophomonas": ("Proteobacteria", "Xanthomonadaceae"),
    "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium": ("Proteobacteria", "Rhizobiaceae"),
}


def natural_taxon_key(x):
    m = re.fullmatch(r"species(\d+)", str(x))
    if m:
        return (0, int(m.group(1)))
    return (1, str(x))


def parse_taxon_levels(taxon):
    levels = {}

    for part in str(taxon or "").split(";"):
        part = part.strip()

        if "__" not in part:
            continue

        prefix, value = part.split("__", 1)
        value = value.strip()

        if value:
            levels[prefix] = value

    return levels


def _find_annotation_file():
    for p in ANNOTATION_CANDIDATES:
        if os.path.exists(p):
            return p
    return None


def _read_annotation_table(path):
    ext = os.path.splitext(path)[1].lower()

    if ext in (".xlsx", ".xlsm"):
        try:
            df = pd.read_excel(path, sheet_name=ANNOTATION_XLSX_SHEET, dtype=str)
        except ValueError:
            df = pd.read_excel(path, sheet_name=0, dtype=str)
    elif ext in (".tsv", ".txt"):
        df = pd.read_csv(path, sep="\t", dtype=str)
    elif ext == ".csv":
        df = pd.read_csv(path, dtype=str)
    else:
        raise ValueError(
            f"Unsupported annotation file extension: {ext!r} ({path}). "
            "Expected .xlsx / .tsv / .csv."
        )

    return df.fillna("")


def load_annotation_map():
    path = _find_annotation_file()

    if path is None:
        tried = "\n    ".join(ANNOTATION_CANDIDATES)
        raise FileNotFoundError(
            "Annotation map not found. Tried:\n    " + tried +
            "\nPut Unified_species_annotations.xlsx or mapping .xlsx/.tsv/.csv under "
            f"{DATA_DIR}/ before running this script."
        )

    ann = _read_annotation_table(path)

    ann = ann.rename(columns={
        "Experiment": "experiment",
        "Species": "taxon",
        "Unified_annotation": "unified_annotation",
        "Genus": "genus",
        "Family": "family",
        "Species_name": "species_name",
        "Original_Taxon": "original_taxon",
    })

    for col in [
        "experiment", "taxon", "unified_annotation",
        "genus", "family", "species_name", "original_taxon"
    ]:
        if col not in ann.columns:
            ann[col] = ""

    ann["experiment"] = ann["experiment"].str.lower()

    parsed = ann["original_taxon"].apply(parse_taxon_levels)
    ann["phylum"] = parsed.apply(lambda x: x.get("p", ""))
    ann["family_from_taxon"] = parsed.apply(lambda x: x.get("f", ""))

    ann.loc[ann["family"] == "", "family"] = ann.loc[
        ann["family"] == "", "family_from_taxon"
    ]

    for idx, row in ann.iterrows():
        genus = row["genus"]

        if genus in GENUS_TAXONOMY:
            phylum, family = GENUS_TAXONOMY[genus]

            if not row["phylum"]:
                ann.at[idx, "phylum"] = phylum

            if not row["family"]:
                ann.at[idx, "family"] = family

    ann.loc[ann["unified_annotation"] == "", "unified_annotation"] = ann["taxon"]
    ann.loc[ann["genus"] == "", "genus"] = ann["unified_annotation"]
    ann.loc[ann["family"] == "", "family"] = ann["genus"]
    ann.loc[ann["phylum"] == "", "phylum"] = "Other"

    keep = [
        "experiment", "taxon", "unified_annotation",
        "genus", "family", "phylum", "species_name"
    ]

    ann = ann[keep].drop_duplicates(["experiment", "taxon"])

    n_other = int((ann["phylum"] == "Other").sum())

    print(f"  -> Loaded annotation map: {path} ({len(ann)} rows)")

    if n_other:
        print(
            f"     [warn] {n_other} rows still fell back to phylum='Other' "
            "after GENUS_TAXONOMY lookup; check their genus names."
        )

    return ann


def expand_annotation_aliases(annotation_df):
    if annotation_df.empty:
        return annotation_df

    alias = annotation_df.copy()
    alias["taxon"] = alias["unified_annotation"]
    alias = alias[alias["taxon"] != ""]
    return (
        pd.concat([annotation_df, alias], ignore_index=True)
        .drop_duplicates(["experiment", "taxon"], keep="first")
    )


def attach_annotations(df, annotation_df):
    cols = ["unified_annotation", "genus", "family", "phylum", "species_name"]

    if annotation_df.empty:
        for col in cols:
            df[col] = ""
        return df

    return df.merge(
        expand_annotation_aliases(annotation_df),
        on=["experiment", "taxon"],
        how="left",
    )


# =====================================================================
#  Color palette
#  Experiment-specific genus-first academic palette
# =====================================================================
#
# 设计哲学
# --------
# 你最初版本能产出"颜色丰富、相邻 species 对比强"的关键在于：
#   1. 平铺所有 species 到完整可用 hue 空间
#   2. interleaving 让相邻位置 hue 距离最大化
#   3. 6 种 sat/lum 组合循环再强化区分度
#
# 之前我做的 phylum/family hue band 分配，会把整个 hue 空间切碎，
# 导致单个 family 内的 species 共享狭窄色相区，颜色趋同。
#
# 这版只做一件事：把 species 在"虚拟 hue 轴"上的排序，
# 从原来的字典序改成分类层级排序 (phylum → family → genus → number)。
# 其他一切（色相空间、interleaving、sat/lum 循环）完全沿用原版。
#
# 结果：
#   - 同 family/genus 的 species 在 hue 轴上位置接近 → 色相相近
#   - 但 interleaving 仍然把它们的实际位置打散 → 仍可清晰区分
#   - 整体颜色丰富度与原版一致

# 可用 hue 区间：跳过红色 + 跳过纯绿色（沿用原版）
HUE_SEGMENTS = [(0.13, 0.25), (0.46, 0.88)]

# Sat/Lum 6 种组合循环（完全沿用原版数值）
GENUS_VARIANTS = [
    (0.000, 0.68, 0.46),
    (0.012, 0.55, 0.58),
    (-0.012, 0.78, 0.40),
    (0.024, 0.50, 0.50),
    (-0.024, 0.72, 0.62),
    (0.036, 0.60, 0.38),
    (-0.036, 0.48, 0.54),
    (0.048, 0.74, 0.48),
]

# mortality / temperature 在虚拟 hue 轴上的起始偏移（沿用原版思路）
HUE_OFFSET = {
    "mortality":   0.0,
    "temperature": 0.37,
}


# Article-inspired taxonomic palette.
# Broad groups follow the visual grammar in the reference trees:
# Proteobacteria: warm oranges/pinks/purples plus a few cyan accents;
# Firmicutes: greens with purple accents for Staphylococcaceae;
# Bacteroidota: blues; Actinobacteriota: grey-blue/teal.
# Species within the same genus are generated as nearby variants of the
# genus base color, so abundance plots and iTOL species strips remain matched.
FAMILY_BASE_COLORS = {
    # Proteobacteria
    "Enterobacteriaceae": "#F28E2B",
    "Erwiniaceae": "#E15759",
    "Burkholderiaceae": "#CC79A7",
    "Pseudomonadaceae": "#00A6A6",
    "Xanthomonadaceae": "#FF8C69",
    "Rhizobiaceae": "#B07AA1",
    "Aeromonadaceae": "#EDC948",
    # Firmicutes
    "Bacillaceae": "#2CA25F",
    "Planococcaceae": "#66BD63",
    "Paenibacillaceae": "#A6D854",
    "Staphylococcaceae": "#7B3294",
    "Exiguobacteraceae": "#1B9E77",
    # Bacteroidota
    "Flavobacteriaceae": "#2F80ED",
    "Weeksellaceae": "#4EA3F1",
    "Sphingobacteriaceae": "#1F5FBF",
    # Actinobacteriota
    "Microbacteriaceae": "#4C7890",
}

PHYLUM_BASE_COLORS = {
    "Proteobacteria": "#F28E2B",
    "Firmicutes": "#2CA25F",
    "Bacteroidota": "#2F80ED",
    "Actinobacteriota": "#4C7890",
    "Other": "#8F8F8F",
}

GENUS_BASE_VARIANTS = [
    (0.000, 1.00, 1.00),
    (0.080, 0.95, 1.16),
    (-0.080, 1.08, 0.86),
    (0.150, 0.90, 1.08),
    (-0.150, 1.12, 0.78),
    (0.220, 1.00, 0.96),
]

SPECIES_VARIANTS = [
    (0.000, 1.00, 1.00),
    (0.018, 0.88, 1.18),
    (-0.018, 1.12, 0.82),
    (0.036, 0.95, 0.96),
    (-0.036, 0.82, 1.30),
    (0.054, 1.15, 0.72),
]


def _hue_total_len():
    return sum(hi - lo for lo, hi in HUE_SEGMENTS)


def _pos_to_hue(pos):
    """虚拟位置 → 真实 hue（跳过红/绿禁区）"""
    total_len = _hue_total_len()
    pos = pos % total_len

    acc = 0.0
    for lo, hi in HUE_SEGMENTS:
        seg_len = hi - lo
        if pos < acc + seg_len:
            return lo + (pos - acc)
        acc += seg_len

    return HUE_SEGMENTS[-1][1]


def _hsl_to_hex(h, s, l):
    h = h % 1.0
    r, g, b = colorsys.hls_to_rgb(h, l, s)
    return "#{:02X}{:02X}{:02X}".format(
        int(round(r * 255)),
        int(round(g * 255)),
        int(round(b * 255)),
    )


def _hex_to_hls(hex_color):
    text = str(hex_color).lstrip("#")
    r = int(text[0:2], 16) / 255.0
    g = int(text[2:4], 16) / 255.0
    b = int(text[4:6], 16) / 255.0
    return colorsys.rgb_to_hls(r, g, b)


def _variant_hex(base_hex, dh=0.0, sat_mult=1.0, light_mult=1.0):
    h, l, s = _hex_to_hls(base_hex)
    s = min(max(s * sat_mult, 0.34), 0.88)
    l = min(max(l * light_mult, 0.28), 0.72)
    return _hsl_to_hex(h + dh, s, l)


def _shift_hue_inside_segments(h, dh):
    h = h % 1.0
    for lo, hi in HUE_SEGMENTS:
        if lo <= h <= hi:
            return min(max(h + dh, lo + 0.006), hi - 0.006)
    return _pos_to_hue(dh)


def _interleave_indices(n):
    """原版 interleaving: 0, n//2, 1, n//2+1, ..."""
    if n <= 1:
        return list(range(n))

    half = n // 2
    order = []

    for i in range(half):
        order.append(i)
        if i + half < n:
            order.append(i + half)

    if len(order) < n:
        for i in range(n):
            if i not in order:
                order.append(i)

    return order


def _generate_palette_for_experiment(df_exp, hue_offset):
    """
    一个实验内的 species 染色。

    文章图式规则：
      1. family/phylum 决定一个离散 base color，而不是连续 hue 轴；
      2. 同 family 内不同 genus 用同一 base 的明显变体；
      3. 同 genus 内不同 species 用更近的明暗/饱和度变体；
      4. mortality / temperature 仍然分开生成，保证两个 library 独立映射。
    """
    df = df_exp.copy()

    # 填默认值，保证排序稳定
    for col in ["phylum", "family", "genus", "unified_annotation"]:
        if col not in df.columns:
            df[col] = ""
        df[col] = df[col].fillna("")

    df.loc[df["unified_annotation"] == "", "unified_annotation"] = df["taxon"]
    df.loc[df["genus"] == "", "genus"] = df["unified_annotation"]
    df.loc[df["family"] == "", "family"] = df["genus"]
    df.loc[df["phylum"] == "", "phylum"] = "Other"

    df["taxon_number"] = df["taxon"].apply(lambda x: natural_taxon_key(x)[1])

    # 关键：按分类层级排序，而不是字典序
    df = df.sort_values(
        ["phylum", "family", "genus", "taxon_number", "taxon"]
    ).reset_index(drop=True)

    if len(df) == 0:
        df["color"] = []
        return df

    group_order = (
        df[["phylum", "family", "genus"]]
        .drop_duplicates()
        .sort_values(["phylum", "family", "genus"])
        .reset_index(drop=True)
    )
    group_order["family_key"] = (
        group_order["phylum"].astype(str) + "||" + group_order["family"].astype(str)
    )
    group_order["genus_rank_in_family"] = group_order.groupby("family_key").cumcount()

    genus_base = {}
    for _, row in group_order.iterrows():
        family = row["family"]
        phylum = row["phylum"]
        base = FAMILY_BASE_COLORS.get(
            family, PHYLUM_BASE_COLORS.get(phylum, PHYLUM_BASE_COLORS["Other"])
        )
        dh, sm, lm = GENUS_BASE_VARIANTS[
            int(row["genus_rank_in_family"]) % len(GENUS_BASE_VARIANTS)
        ]
        genus_base[row["genus"]] = _variant_hex(base, dh, sm, lm)

    species_rank = df.groupby("genus").cumcount()

    colors = []
    for i, row in df.iterrows():
        dh, sm, lm = SPECIES_VARIANTS[
            int(species_rank.loc[i]) % len(SPECIES_VARIANTS)
        ]
        colors.append(_variant_hex(genus_base[row["genus"]], dh, sm, lm))

    df["color"] = colors
    return df


def generate_color_palette(mortality_taxa, temperature_taxa, annotation_df=None):
    """
    mortality / temperature 独立染色，颜色不共享。
    返回长格式 color_map: taxon, experiment, color, + annotation 字段
    """
    rows = []

    for taxon in mortality_taxa:
        rows.append({"taxon": taxon, "experiment": "mortality"})
    for taxon in temperature_taxa:
        rows.append({"taxon": taxon, "experiment": "temperature"})

    df = pd.DataFrame(rows)

    if annotation_df is not None and not annotation_df.empty:
        df = df.merge(
            expand_annotation_aliases(annotation_df),
            on=["experiment", "taxon"],
            how="left",
        )

    for col in ["unified_annotation", "genus", "family", "phylum", "species_name"]:
        if col not in df.columns:
            df[col] = ""
        df[col] = df[col].fillna("")

    df.loc[df["unified_annotation"] == "", "unified_annotation"] = df["taxon"]
    df.loc[df["genus"] == "", "genus"] = df["unified_annotation"]
    df.loc[df["family"] == "", "family"] = df["genus"]
    df.loc[df["phylum"] == "", "phylum"] = "Other"

    out = []
    for exp in ["mortality", "temperature"]:
        df_exp = df[df["experiment"] == exp].copy()
        if df_exp.empty:
            continue
        out.append(
            _generate_palette_for_experiment(df_exp, hue_offset=HUE_OFFSET[exp])
        )

    color_df = pd.concat(out, ignore_index=True)

    color_df["display_label"] = (
        color_df["taxon"] + " (" + color_df["unified_annotation"] + ")"
    )
    color_df["color_group"] = color_df["phylum"] + " | " + color_df["family"]

    color_df["taxon_number"] = color_df["taxon"].apply(
        lambda x: natural_taxon_key(x)[1]
    )

    color_df = (
        color_df
        .sort_values(["experiment", "phylum", "family", "genus", "taxon_number", "taxon"])
        .drop(columns=["taxon_number"])
        .reset_index(drop=True)
    )

    keep = [
        "taxon", "experiment", "color", "display_label",
        "unified_annotation", "genus", "family", "phylum",
        "species_name", "color_group",
    ]
    return color_df[keep]


# =====================================================================
#  Main
# =====================================================================

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    all_dfs = []
    annotation_df = load_annotation_map()

    for exp in ["mortality", "temperature"]:
        print(f"\n处理 {exp} (OD背景: {OD_BACKGROUND[exp]})...")

        seq_df = parse_seq_excel(SEQ_FILES[exp], exp)
        od_df  = parse_od_excel(OD_FILES[exp], exp)

        df = compute_abundance(seq_df, od_df)
        df = attach_annotations(df, annotation_df)

        print(f"  -> 记录数: {len(df)}, Taxon: {df['taxon'].nunique()}")

        for cond in sorted(df["condition"].unique()):
            cond_dir = os.path.join(OUT_DIR, exp, cond)
            os.makedirs(cond_dir, exist_ok=True)

            out_path = os.path.join(cond_dir, "abs_abundance.csv")

            (
                df[df["condition"] == cond]
                .sort_values(["community", "taxon", "replica", "day"])
                .to_csv(out_path, index=False)
            )

        print(f"  -> 已保存至 {OUT_DIR}/{exp}/W1~W5/")

        all_dfs.append(df)

    print("\n生成色板（两个 library 独立，genus-first academic palette）...")

    all_df = pd.concat(all_dfs, ignore_index=True)

    m_taxa = set(all_df[all_df["experiment"] == "mortality"]["taxon"])
    t_taxa = set(all_df[all_df["experiment"] == "temperature"]["taxon"])

    color_df = generate_color_palette(m_taxa, t_taxa, annotation_df)

    color_path = os.path.join(OUT_DIR, "color_map.csv")
    color_df.to_csv(color_path, index=False)

    print(f"  -> {len(color_df)} 色 → {color_path}")

    print("\n完成！")


if __name__ == "__main__":
    main()
