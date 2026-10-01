"""
Plot total biomass fluctuation against relative-abundance community CV.

X axis:
    total_biomass_std = SD(OD_t) within the analysis window

Y axis:
    community_cv_rel = sum_i SD(relative_abundance_i_clean) / 100

The script first looks for community_cv_rel in the formal fluctuation output.
If that column is not present yet, it falls back to the temporary diagnostic
table produced during method exploration.
"""

from pathlib import Path
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


BASE_DIR = Path(os.environ.get("COMMUNITY_WORKSPACE", Path(__file__).resolve().parent))
if os.name == "nt" or not (BASE_DIR / "processed").exists():
    BASE_DIR = Path(__file__).resolve().parent

PROC_DIR = BASE_DIR / "processed"
FIG_DIR = BASE_DIR / "figures" / "total_biomass_std_vs_rel_cv"
OUT_DIR = PROC_DIR / "total_biomass_std_vs_rel_cv"
FIG_DIR.mkdir(parents=True, exist_ok=True)
OUT_DIR.mkdir(parents=True, exist_ok=True)

MAIN_WINDOWS = {
    "temperature": "full",
    "mortality": "early",
}

CV_THRESHOLD = 0.265
BC_MAX_THRESHOLD = 0.20

COLORS = {
    "stable": "#9B8EC4",
    "cv_high": "#F4A460",
    "bcmax_high": "#4C78A8",
    "both_high": "#55A868",
}


def read_main_windows() -> pd.DataFrame:
    frames = []
    for experiment, window in MAIN_WINDOWS.items():
        path = PROC_DIR / "fluctuations" / window / "community_level.csv"
        if not path.exists():
            continue
        df = pd.read_csv(path)
        df["window"] = window
        frames.append(df[df["experiment"] == experiment].copy())
    if not frames:
        raise FileNotFoundError("No main fluctuation community_level.csv files found.")
    return pd.concat(frames, ignore_index=True)


def attach_relative_cv(df: pd.DataFrame) -> pd.DataFrame:
    if "community_cv_rel" in df.columns:
        return df

    fallback = PROC_DIR / "statistical_tests" / "community_cv_relative_abundance_windows.csv"
    if not fallback.exists():
        raise FileNotFoundError(
            "community_cv_rel is absent from community_level.csv and fallback "
            f"table was not found: {fallback}"
        )

    rel = pd.read_csv(fallback)
    keep_cols = [
        "window",
        "experiment",
        "condition",
        "community",
        "replica",
        "community_cv_rel",
        "sum_rel_std",
        "temporal_bc_max_calc",
    ]
    keep_cols = [c for c in keep_cols if c in rel.columns]
    return df.merge(
        rel[keep_cols],
        on=["window", "experiment", "condition", "community", "replica"],
        how="left",
    )


def prepare_data() -> pd.DataFrame:
    df = attach_relative_cv(read_main_windows())
    if "temporal_bc_max" not in df.columns and "temporal_bc_max_calc" in df.columns:
        df["temporal_bc_max"] = df["temporal_bc_max_calc"]

    numeric_cols = [
        "community_cv",
        "community_cv_rel",
        "total_biomass_std",
        "temporal_bc_max",
    ]
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    if "collapsed" not in df.columns:
        df["collapsed"] = False
    df["collapsed"] = df["collapsed"].fillna(False).astype(bool)

    df = df[
        np.isfinite(df["community_cv_rel"])
        & np.isfinite(df["total_biomass_std"])
    ].copy()
    df["cv_high"] = (df["community_cv"] >= CV_THRESHOLD) & (~df["collapsed"])
    df["bcmax_high"] = (
        (df.get("temporal_bc_max", np.nan) >= BC_MAX_THRESHOLD)
        & (~df["collapsed"])
    )
    df["instability_class"] = np.select(
        [
            df["cv_high"] & df["bcmax_high"],
            df["cv_high"],
            df["bcmax_high"],
        ],
        ["both_high", "cv_high", "bcmax_high"],
        default="stable",
    )
    return df


def save_input_table(df: pd.DataFrame) -> None:
    cols = [
        "experiment",
        "window",
        "condition",
        "community",
        "replica",
        "total_biomass_std",
        "community_cv_rel",
        "community_cv",
        "temporal_bc_max",
        "collapsed",
        "instability_class",
    ]
    cols = [c for c in cols if c in df.columns]
    df[cols].to_csv(OUT_DIR / "total_biomass_std_vs_community_cv_rel_values.csv",
                    index=False)


def scatter_panel(ax, df: pd.DataFrame, title: str) -> None:
    for cls, color in COLORS.items():
        sub = df[df["instability_class"] == cls]
        if sub.empty:
            continue
        ax.scatter(
            sub["total_biomass_std"],
            sub["community_cv_rel"],
            s=22,
            alpha=0.82,
            c=color,
            edgecolors="none",
            label=cls,
        )
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("Total biomass std", fontsize=9)
    ax.set_ylabel("Relative-abundance community CV", fontsize=9)
    ax.tick_params(labelsize=8, direction="in", top=True, right=True)
    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(True)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)


def main() -> None:
    df = prepare_data()
    save_input_table(df)

    for experiment, sub in df.groupby("experiment"):
        fig, ax = plt.subplots(figsize=(3.2, 3.0), dpi=300)
        scatter_panel(ax, sub, experiment)
        ax.legend(frameon=False, fontsize=7, loc="best")
        fig.tight_layout()
        fig.savefig(FIG_DIR / f"{experiment}_total_biomass_std_vs_community_cv_rel.pdf")
        fig.savefig(FIG_DIR / f"{experiment}_total_biomass_std_vs_community_cv_rel.png")
        plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(6.4, 3.0), dpi=300, sharey=False)
    for ax, experiment in zip(axes, ["mortality", "temperature"]):
        sub = df[df["experiment"] == experiment]
        scatter_panel(ax, sub, experiment)
    handles, labels = axes[-1].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, fontsize=7,
               loc="center right", bbox_to_anchor=(1.02, 0.5))
    fig.tight_layout(rect=(0, 0, 0.88, 1))
    fig.savefig(FIG_DIR / "main_total_biomass_std_vs_community_cv_rel_facets.pdf")
    fig.savefig(FIG_DIR / "main_total_biomass_std_vs_community_cv_rel_facets.png")
    plt.close(fig)

    print(f"Saved plots to {FIG_DIR}")
    print(f"Saved values to {OUT_DIR}")


if __name__ == "__main__":
    main()
