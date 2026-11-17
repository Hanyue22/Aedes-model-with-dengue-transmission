#!/usr/bin/env python3
"""
Compare hl vs hs simulation results.

Single combined figure:
  - Top row:  Drive Homozygotes, DR, R2 (time-series)
  - Bottom row: Total New Infections (scatter + mean ± SEM)
  - Single legend at bottom.

Filename format:  mod_{type}_rep{N}_{timestamps}.xlsx

Usage:
    python plot_hl_vs_hs.py
"""

import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

# ── Configuration ──────────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).resolve().parent

# Which time-series metrics to plot (top row)
METRICS_TS = [
    ("Rate_homo_fe", "Drive Homozygotes"),
    ("DR_1",         "Drive Allele"),
    ("R2_1",         "R2 Allele"),
]

START_WEEK = 16
YEARS = 20
END_WEEK = START_WEEK + YEARS * 52 - 1

# Legend display names
LEGEND_LABEL = {
    "hl": "Haplolethal target",
    "hs": "Haplosufficient target",
}

# ── Visual parameters ──────────────────────────────────────────────────────
COLOR = {"hl": "#2166ac", "hs": "#d6604d"}
DPI = 300

# Font sizes
TITLE_FONTSIZE = 16
SUBTITLE_FONTSIZE = 16
AXIS_LABEL_FONTSIZE = 16
TICK_LABELSIZE = 14
LEGEND_FONTSIZE = 14

# Time-series subplots
TS_SCATTER_ALPHA = 0.18
TS_SCATTER_SIZE = 3
TS_LINE_WIDTH = 1.6
TS_SEM_ALPHA = 0.2

# Total-infections subplot
TI_SCATTER_SIZE = 50
TI_SCATTER_ALPHA = 0.7
TI_LINE_WIDTH = 2.5
TI_CAPSIZE = 10
TI_JITTER = 0.2
TI_BAR_HW = 0.28

# Figure layout
FIG_WIDTH = 10
FIG_HEIGHT = 10
HEIGHT_RATIOS = [1, 1.8]   # bottom panel 1.5× taller
HSPACE = 0.5              # vertical gap between rows
WSPACE = 0.35              # horizontal gap between columns


# ── File discovery ─────────────────────────────────────────────────────────
def parse_filename(filename: str) -> dict | None:
    m = re.match(r"mod_(hl|hs)_rep(\d+)_.*\.xlsx", filename)
    if not m:
        return None
    return {"type": m.group(1), "rep": int(m.group(2))}


def discover_files(directory: Path) -> pd.DataFrame:
    records = []
    for f in sorted(directory.glob("mod_*.xlsx")):
        info = parse_filename(f.name)
        if info is None:
            print(f"  ⚠ skipped: {f.name}")
            continue
        info["path"] = f
        records.append(info)

    df = pd.DataFrame(records)
    for typ in ["hl", "hs"]:
        sub = df[df["type"] == typ]
        print(f"  {typ}: {len(sub)} files, {sub['rep'].nunique()} unique reps")
    return df


# ── Data loading ───────────────────────────────────────────────────────────
def load_data(file_df: pd.DataFrame) -> list[dict]:
    col_map = {col: label for col, label in METRICS_TS}
    col_map["Native New Infections"] = "Weekly New Infections"
    all_cols = (["Generation", "Native New Infections"]
                + [c for c, _ in METRICS_TS])

    records = []
    for _, row in file_df.iterrows():
        df = pd.read_excel(row["path"], header=0, usecols=all_cols)
        df = df.set_index("Generation")
        df = df.rename(columns=col_map)
        records.append({"type": row["type"], "rep": row["rep"], "df": df})
    print(f"  Loaded {len(records)} files")
    return records


# ── Helpers ────────────────────────────────────────────────────────────────
def align_matrix(dfs: list[pd.DataFrame], metric: str
                 ) -> tuple[np.ndarray, np.ndarray]:
    all_gens = sorted(set().union(*(df.index for df in dfs)))
    rows = []
    for df in dfs:
        s = df.get(metric, None)
        if s is None:
            continue
        rows.append(s.reindex(all_gens).values)
    return np.array(all_gens), np.array(rows)


def mean_sem(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    n = matrix.shape[0]
    mean = np.nanmean(matrix, axis=0)
    sem = np.nanstd(matrix, axis=0, ddof=1) / np.sqrt(n)
    return mean, sem


# ── Combined figure ────────────────────────────────────────────────────────
def plot_combined(data: list[dict], output_path: Path):
    """Top: 3 time-series panels. Bottom: Total New Infections.
    Single legend at the bottom."""

    fig = plt.figure(figsize=(FIG_WIDTH, FIG_HEIGHT))
    gs = fig.add_gridspec(
        2, 3,
        height_ratios=HEIGHT_RATIOS,
        hspace=HSPACE, wspace=WSPACE,
        left=0.08, right=0.97, top=0.95, bottom=0.10,
    )

    # ── Row 0: time-series ────────────────────────────────────────────────
    ax_left = None   # leftmost top-row axes, for y-label alignment
    ax_legend = None  # center top-row axes, for legend anchor
    for idx, (col_name, label) in enumerate(METRICS_TS):
        ax = fig.add_subplot(gs[0, idx])
        if idx == 0:
            ax_left = ax
        if idx == 1:
            ax_legend = ax
        is_homo = (col_name == "Rate_homo_fe")
        for typ in ["hl", "hs"]:
            dfs = [r["df"] for r in data if r["type"] == typ]
            if not dfs:
                continue
            gens, mat = align_matrix(dfs, label)
            if mat.size == 0:
                continue

            if is_homo:
                mat = mat / 100.0

            for row_idx in range(mat.shape[0]):
                ax.scatter(
                    gens, mat[row_idx],
                    s=TS_SCATTER_SIZE, color=COLOR[typ],
                    alpha=TS_SCATTER_ALPHA, linewidths=0, zorder=1,
                )
            mu, sem = mean_sem(mat)
            ax.fill_between(
                gens, mu - sem, mu + sem,
                color=COLOR[typ], alpha=TS_SEM_ALPHA, linewidth=0, zorder=2,
            )
            ax.plot(gens, mu, color=COLOR[typ], linewidth=TS_LINE_WIDTH,
                    zorder=3)

        ax.set_title(label, fontsize=SUBTITLE_FONTSIZE, fontweight="bold")
        ax.set_xlabel("Week", fontsize=AXIS_LABEL_FONTSIZE, fontweight="bold")
        ax.set_ylim(0, 1.05)
        ax.tick_params(labelsize=TICK_LABELSIZE)

        if is_homo:
            ax.yaxis.set_major_formatter(
                mticker.FuncFormatter(lambda v, _: f"{v:.1f}"))
            ax.set_ylabel("Frequency", fontsize=AXIS_LABEL_FONTSIZE, fontweight="bold")
        else:
            ax.yaxis.set_major_formatter(
                mticker.FuncFormatter(lambda v, _: f"{v:.1f}")
            )

    # ── Row 1: Total New Infections (span all 3 cols) ─────────────────────
    ax_ti = fig.add_subplot(gs[1, :])

    metric = "Weekly New Infections"
    rng = np.random.default_rng(42)
    x_positions = {"hl": 0, "hs": 1}

    for typ, x in x_positions.items():
        totals = []
        for rec in data:
            if rec["type"] != typ:
                continue
            subset = rec["df"].loc[START_WEEK:END_WEEK, metric]
            totals.append(subset.sum())

        totals_arr = np.array(totals)
        if len(totals_arr) == 0:
            continue

        mean_val = np.mean(totals_arr)
        sem_val = np.std(totals_arr, ddof=1) / np.sqrt(len(totals_arr))

        jitter = rng.uniform(-TI_JITTER, TI_JITTER, size=len(totals_arr))
        ax_ti.scatter(
            x + jitter, totals_arr,
            s=TI_SCATTER_SIZE, color=COLOR[typ], alpha=TI_SCATTER_ALPHA,
            linewidths=0.5, edgecolors="white", zorder=3,
        )

        ax_ti.plot(
            [x - TI_BAR_HW, x + TI_BAR_HW], [mean_val, mean_val],
            color=COLOR[typ], linewidth=TI_LINE_WIDTH, linestyle="-", zorder=4,
        )
        ax_ti.errorbar(
            x, mean_val, yerr=sem_val, color=COLOR[typ],
            linewidth=TI_LINE_WIDTH, capsize=TI_CAPSIZE, capthick=2, zorder=4,
        )

    ax_ti.set_xticks([0, 1])
    ax_ti.set_xticklabels(
        ["Haplolethal target", "Haplosufficient target"],
        fontsize=SUBTITLE_FONTSIZE, fontweight="bold",
    )
    ax_ti.set_xlim(-0.5, 1.5)
    ax_ti.set_ylim(0, 1400)
    ax_ti.set_ylabel("Total New Infections", fontsize=AXIS_LABEL_FONTSIZE, fontweight="bold")
    ax_ti.set_title(
        f"Total New Infections (20 Years)",
        fontsize=SUBTITLE_FONTSIZE, fontweight="bold",
    )
    ax_ti.tick_params(labelsize=TICK_LABELSIZE)
    ax_ti.yaxis.set_major_formatter(
        mticker.FuncFormatter(lambda v, _: f"{v:.0f}")
    )

    # ── Align y-axis labels (same X, each axis's own Y center) ──────────
    from matplotlib.transforms import blended_transform_factory
    fig.canvas.draw()
    label_x = ax_left.get_position().x0 - 0.06
    trans_left = blended_transform_factory(fig.transFigure, ax_left.transAxes)
    trans_ti = blended_transform_factory(fig.transFigure, ax_ti.transAxes)
    ax_left.yaxis.set_label_coords(label_x, 0.5, transform=trans_left)
    ax_ti.yaxis.set_label_coords(label_x, 0.5, transform=trans_ti)

    # ── Unified legend between top and bottom rows ───────────────────────
    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], color=COLOR["hl"], linewidth=2.5,
               label=LEGEND_LABEL["hl"]),
        Line2D([0], [0], color=COLOR["hs"], linewidth=2.5,
               label=LEGEND_LABEL["hs"]),
    ]
    ax_legend.legend(
        handles=handles, loc="upper center",
        bbox_to_anchor=(0.5, -0.22),
        fontsize=LEGEND_FONTSIZE, ncol=2, frameon=False,
    )

    fig.savefig(output_path, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  ✓ Saved {output_path}")


# ── Main ───────────────────────────────────────────────────────────────────
def main():
    print("=" * 60)
    print("hl vs hs comparison")
    print(f"Working directory: {SCRIPT_DIR}")
    print(f"Window: week {START_WEEK}–{END_WEEK} ({YEARS} yrs)")
    print("=" * 60)

    print("\n[1/3] Scanning for xlsx files …")
    file_df = discover_files(SCRIPT_DIR)
    if file_df.empty:
        print("No files found. Exiting.")
        return

    print("\n[2/3] Loading data …")
    data = load_data(file_df)

    print("\n[3/3] Plotting …")
    plot_combined(data, SCRIPT_DIR / "combined.png")

    print("\n✅ Done.")


if __name__ == "__main__":
    main()
