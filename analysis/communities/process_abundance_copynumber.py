"""
process_abundance.py
主目录: /home/hachi/Tem_mortality_workspace

功能
----
1. 从统一测序 Excel 和 OD Excel 中解析数据
2. 加载 species annotation，并读取 16S rRNA gene copy number
3. 基于 16S rRNA gene copy number 校正 reads
4. 计算 copy-number corrected relative abundance 和 absolute abundance
5. 同时保留 raw relative abundance / raw absolute abundance 作为对照
6. 输出 processed/{experiment}/W*/abs_abundance.csv
7. 输出 processed/color_map.csv

重要说明
--------
请把带 16S copy number 的注释表放在：

    /home/hachi/Tem_mortality_workspace/data/Unified_species_annotations_with_16S_copy_number.xlsx

脚本会优先读取这个文件。
如果没有找到它，才会依次尝试读取原来的 Unified_species_annotations.xlsx 等文件。

16S 多拷贝校正逻辑
------------------
原始相对丰度：

    rel_abund_raw = reads / sum(reads) * 100

多拷贝校正：

    reads_copy_corrected = reads / copy_number_16s
    rel_abund = reads_copy_corrected / sum(reads_copy_corrected) * 100
    abs_abund = rel_abund / 100 * OD

也就是说，输出中的：

    rel_abund_raw  = 未做多拷贝校正的相对丰度
    abs_abund_raw  = 未做多拷贝校正的绝对丰度
    rel_abund      = 做过 16S 多拷贝校正后的相对丰度
    abs_abund      = 做过 16S 多拷贝校正后的绝对丰度

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

# 优先读取带 16S copy number 的注释文件
ANNOTATION_CANDIDATES = [
    os.path.join(DATA_DIR, "Unified_species_annotations_with_16S_copy_number.xlsx"),
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
    """
    计算丰度。

    这里假设 seq_df 已经通过 attach_annotations() 合并了 annotation，
    因而包含 copy_number_16s 字段。

    输出字段说明
    ----------
    reads:
        原始测序 reads

    total_reads:
        每个 sample 的原始 reads 总和

    rel_abund_raw:
        未做 16S 多拷贝校正的相对丰度，单位 %

    copy_number_16s:
        每个 taxon 的 16S rRNA gene copy number

    reads_copy_corrected:
        reads / copy_number_16s

    total_reads_copy_corrected:
        每个 sample 内 copy-number corrected reads 的总和

    rel_abund:
        做过 16S 多拷贝校正后的相对丰度，单位 %

    abs_abund_raw:
        未做 16S 多拷贝校正的绝对丰度，计算为 rel_abund_raw / 100 * OD

    abs_abund:
        做过 16S 多拷贝校正后的绝对丰度，计算为 rel_abund / 100 * OD
    """
    df = seq_df.copy()

    if df.empty:
        return df

    # ------------------------------------------------------------
    # 1. 原始 reads 总和，用于保留 raw relative abundance
    # ------------------------------------------------------------
    raw_totals = (
        df.dropna(subset=["reads"])
        .groupby(["experiment", "condition", "community", "replica", "day"])["reads"]
        .sum()
        .reset_index()
        .rename(columns={"reads": "total_reads"})
    )

    df = df.merge(
        raw_totals,
        on=["experiment", "condition", "community", "replica", "day"],
        how="left",
    )

    df["rel_abund_raw"] = np.where(
        df["reads"].isna() | (df["total_reads"] == 0),
        np.nan,
        df["reads"] / df["total_reads"] * 100,
    )

    # ------------------------------------------------------------
    # 2. 读取 16S copy number
    #    如果没有 copy number，则默认 1.0，相当于不校正
    # ------------------------------------------------------------
    if "copy_number_16s" not in df.columns:
        df["copy_number_16s"] = 1.0

    df["copy_number_16s"] = pd.to_numeric(
        df["copy_number_16s"],
        errors="coerce",
    )

    df["copy_number_16s"] = df["copy_number_16s"].replace(
        [np.inf, -np.inf],
        np.nan,
    )

    # copy number 缺失、为 0 或负数时，默认按 1 处理
    df["copy_number_16s"] = np.where(
        df["copy_number_16s"].isna() | (df["copy_number_16s"] <= 0),
        1.0,
        df["copy_number_16s"],
    )

    # ------------------------------------------------------------
    # 3. 多拷贝校正后的 reads
    # ------------------------------------------------------------
    df["reads_copy_corrected"] = np.where(
        df["reads"].isna(),
        np.nan,
        df["reads"] / df["copy_number_16s"],
    )

    corrected_totals = (
        df.dropna(subset=["reads_copy_corrected"])
        .groupby(["experiment", "condition", "community", "replica", "day"])["reads_copy_corrected"]
        .sum()
        .reset_index()
        .rename(columns={"reads_copy_corrected": "total_reads_copy_corrected"})
    )

    df = df.merge(
        corrected_totals,
        on=["experiment", "condition", "community", "replica", "day"],
        how="left",
    )

    # ------------------------------------------------------------
    # 4. 校正后的 relative abundance
    #    注意：这里 rel_abund 是 copy-number corrected 结果
    # ------------------------------------------------------------
    df["rel_abund"] = np.where(
        df["reads_copy_corrected"].isna() | (df["total_reads_copy_corrected"] == 0),
        np.nan,
        df["reads_copy_corrected"] / df["total_reads_copy_corrected"] * 100,
    )

    # ------------------------------------------------------------
    # 5. 合并 OD，计算 absolute abundance
    # ------------------------------------------------------------
    df = df.merge(
        od_df,
        on=["experiment", "condition", "community", "replica", "day"],
        how="left",
    )

    # 未校正 absolute abundance，保留作对照
    df["abs_abund_raw"] = np.where(
        df["rel_abund_raw"].isna() | df["OD"].isna(),
        np.nan,
        df["rel_abund_raw"] / 100.0 * df["OD"],
    )

    # 校正后的 absolute abundance
    df["abs_abund"] = np.where(
        df["rel_abund"].isna() | df["OD"].isna(),
        np.nan,
        df["rel_abund"] / 100.0 * df["OD"],
    )

    # ------------------------------------------------------------
    # 6. round
    # ------------------------------------------------------------
    df["rel_abund_raw"] = df["rel_abund_raw"].round(6)
    df["rel_abund"] = df["rel_abund"].round(6)

    df["reads_copy_corrected"] = df["reads_copy_corrected"].round(8)
    df["total_reads_copy_corrected"] = df["total_reads_copy_corrected"].round(8)

    df["abs_abund_raw"] = df["abs_abund_raw"].round(8)
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
            "\nPut Unified_species_annotations_with_16S_copy_number.xlsx "
            "or Unified_species_annotations.xlsx under "
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

        # 16S copy-number correction fields
        "16S_copy_number": "copy_number_16s",
        "Copy_number_match_level": "copy_number_match_level",
        "Copy_number_source_taxon": "copy_number_source_taxon",
        "Copy_number_confidence": "copy_number_confidence",
        "Copy_number_note": "copy_number_note",
    })

    for col in [
        "experiment", "taxon", "unified_annotation",
        "genus", "family", "species_name", "original_taxon",
        "copy_number_16s", "copy_number_match_level",
        "copy_number_source_taxon", "copy_number_confidence",
        "copy_number_note",
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
        "genus", "family", "phylum", "species_name",
        "copy_number_16s", "copy_number_match_level",
        "copy_number_source_taxon", "copy_number_confidence",
        "copy_number_note",
    ]

    ann = ann[keep].drop_duplicates(["experiment", "taxon"])

    n_other = int((ann["phylum"] == "Other").sum())
    n_missing_cn = int((pd.to_numeric(ann["copy_number_16s"], errors="coerce").isna()).sum())

    print(f"  -> Loaded annotation map: {path} ({len(ann)} rows)")

    if n_other:
        print(
            f"     [warn] {n_other} rows still fell back to phylum='Other' "
            "after GENUS_TAXONOMY lookup; check their genus names."
        )

    if n_missing_cn:
        print(
            f"     [warn] {n_missing_cn} rows have missing/invalid copy_number_16s; "
            "these rows will be treated as copy_number_16s = 1.0."
        )

    return ann


def attach_annotations(df, annotation_df):
    cols = [
        "unified_annotation", "genus", "family", "phylum", "species_name",
        "copy_number_16s", "copy_number_match_level",
        "copy_number_source_taxon", "copy_number_confidence",
        "copy_number_note",
    ]

    if annotation_df.empty:
        for col in cols:
            df[col] = ""
        return df

    out = df.merge(annotation_df, on=["experiment", "taxon"], how="left")

    for col in cols:
        if col not in out.columns:
            out[col] = ""

    return out


# =====================================================================
#  Color palette
#  原版风格 + 轻量分类排序约束
# =====================================================================

# 可用 hue 区间：跳过红色 + 跳过纯绿色（沿用原版）
HUE_SEGMENTS = [(0.06, 0.28), (0.42, 0.92)]

# Sat/Lum 6 种组合循环（完全沿用原版数值）
SAT_LUM = [
    (0.80, 0.48),
    (0.62, 0.58),
    (0.85, 0.46),
    (0.58, 0.56),
    (0.75, 0.44),
    (0.68, 0.60),
]

# mortality / temperature 在虚拟 hue 轴上的起始偏移（沿用原版思路）
HUE_OFFSET = {
    "mortality":   0.0,
    "temperature": 0.5,
}


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
    r, g, b = colorsys.hls_to_rgb(h, l, s)
    return "#{:02X}{:02X}{:02X}".format(
        int(round(r * 255)),
        int(round(g * 255)),
        int(round(b * 255)),
    )


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

    流程：
      1. 按 phylum → family → genus → taxon_number 排序 species
      2. 在虚拟 hue 轴上等距分配 n 个位置
      3. 用原版 interleaving 把位置打散
      4. sat/lum 6 组循环
      5. hue_offset 让 mortality / temperature 在色环上错开
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

    n = len(df)
    if n == 0:
        df["color"] = []
        return df

    total_len = _hue_total_len()

    # 等距虚拟位置（与原版一致）
    positions = [total_len * i / n for i in range(n)]

    # 原版 interleaving：把相邻 hue slot 打散
    order = _interleave_indices(n)
    interleaved = [positions[i] for i in order]

    # 真实 hue + offset
    hues = [_pos_to_hue(p + hue_offset * total_len) for p in interleaved]

    # sat/lum 6 组循环（slot 顺序，与原版一致）
    colors = []
    for i, h in enumerate(hues):
        s, l = SAT_LUM[i % len(SAT_LUM)]
        colors.append(_hsl_to_hex(h, s, l))

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
        df = df.merge(annotation_df, on=["experiment", "taxon"], how="left")

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

        # 关键修改：先合并 species annotation 和 16S copy number
        seq_df = attach_annotations(seq_df, annotation_df)

        # 再计算 copy-number corrected abundance
        df = compute_abundance(seq_df, od_df)

        print(f"  -> 记录数: {len(df)}, Taxon: {df['taxon'].nunique()}")

        # 检查是否还有未匹配 copy number 的 taxon
        missing_cn = (
            df[["taxon", "copy_number_16s"]]
            .drop_duplicates()
            .assign(copy_number_16s_num=lambda x: pd.to_numeric(x["copy_number_16s"], errors="coerce"))
        )
        n_missing_cn = int(missing_cn["copy_number_16s_num"].isna().sum())
        if n_missing_cn:
            print(
                f"     [warn] {n_missing_cn} taxa have missing/invalid copy_number_16s "
                "and were treated as 1.0."
            )

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

    print("\n生成色板（原版 HSL 风格 + 分类排序约束）...")

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
