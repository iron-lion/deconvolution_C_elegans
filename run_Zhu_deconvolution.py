"""
run_zhu_deconvolution.py
========================
Test all deconvolution algorithms on the real Zhu et al. (2024)
single-worm proteomics dataset with multiple parameter configurations.

Reproduces the full pipeline from notebooks/deconvolution_for_pub.ipynb:
  - Same data loading (zhu_ref_data)
  - Same reference-matrix construction and gene-alignment
  - Same tissue-group aggregation via marker_genes_convert
  - Same CSV and SVG output format

Results are written to  tests/results/
with filenames encoding the algorithm and parameter set used.

Usage (from the project root):
    python tests/run_zhu_deconvolution.py
"""

import os
import sys
import warnings

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")                          # headless — no display needed
import matplotlib.pyplot as plt

# ── path setup ──────────────────────────────────────────────────────────────
ROOT     = os.path.dirname((os.path.abspath(__file__)))
SRC_DIR  = os.path.join(ROOT, "src")
DATA_DIR = os.path.join(ROOT, "data")
OUT_DIR  = os.path.join(ROOT, "results", "Zhu")
print(OUT_DIR)
os.makedirs(OUT_DIR, exist_ok=True)

sys.path.insert(0, SRC_DIR)

from utils import zhu_ref_data, load_cell_to_group
from lsvr import run_multiprocess_deconvolution


# ── colour palette (matches notebook) ───────────────────────────────────────
import seaborn as sns

COLOR_TARGET = [
    "Germline", "Glia", "Gonadal", "Intestine", "Muscle",
    "Neurons", "Uterine-vulval cells", "vulval cells and the epidermis",
]

def _build_color_map(target_columns):
    palette  = sns.color_palette("Set1", 8)
    grays    = sns.color_palette("gist_yarg", 15)
    color_map, ci, gi = [], 0, 0
    for col in target_columns:
        if col in COLOR_TARGET:
            color_map.append(palette[ci]); ci += 1
        else:
            color_map.append(grays[(7 * gi) % 10]); gi += 1
    return color_map


# ── core pipeline ────────────────────────────────────────────────────────────

def build_aligned_matrices(target_df, min_markerscore=0.08, min_marker=5):
    """
    Align the bulk proteomics matrix and the reference matrix to their
    shared marker genes.  Mirrors run_one_lsvr() from the notebook exactly.
    """
    markers = pd.read_csv(
        os.path.join(DATA_DIR, "top_100_markers_postsub.csv"),
        index_col=None, header=0,
    )
    markers = markers[markers["marker_score"] > min_markerscore]

    bulk = target_df[[c for c in target_df.columns if len(str(c)) > 0]].copy()
    bulk = bulk.T.reindex(set(markers["gene_id"])).T
    bulk = bulk.dropna(axis=1, how="all")

    reference_df = markers.pivot_table(
        index="gene_id", columns="cell_group", values="mean_expression"
    )
    reference_df = reference_df.reindex(
        bulk.columns.intersection(reference_df.index)
    ).fillna(0)
    reference_df = reference_df.loc[
        :, (reference_df > 0).sum(axis=0) > min_marker
    ]

    # Keep only the genes present in the (filtered) reference
    bulk = bulk[reference_df.index]

    print(f"  Bulk matrix  : {bulk.shape[0]} samples × {bulk.shape[1]} genes")
    print(f"  Reference    : {reference_df.shape[0]} genes × {reference_df.shape[1]} cell types")
    return bulk, reference_df


def aggregate_to_tissues(raw_prop_df, marker_genes_convert):
    """
    Map fine-grained cell types → tissue groups and sum proportions.
    Mirrors the notebook's post-processing step.
    """
    agg = raw_prop_df.T.copy()
    agg["group"] = agg.index.map(marker_genes_convert)
    agg["group"] = [
        str(g).split(" is ")[0].split(" are ")[0]
        if isinstance(g, str) and "These gene" not in g
        else g
        for g in agg["group"]
    ]
    agg = agg.groupby("group").sum()
    return agg.T


def save_results(
    prop_df,           # samples × cell-types  (NO timepoints column)
    time_points,       # list of timepoint values, one per sample
    tag,               # file-name identifier
    marker_genes_convert,
):
    """
    Save the raw proportion CSV and the tissue-aggregated plot-data CSV + SVG.
    timepoints are used only for day-level aggregation and never saved as a
    proportion column or plotted as a cell type.
    """
    # 1. Raw proportions (cell-type level only — no timepoints column)
    raw_csv = os.path.join(OUT_DIR, f"Zhu_{tag}_proportion.csv")
    prop_df.to_csv(raw_csv)
    print(f"    Saved: {raw_csv}")

    # 2. Aggregate to tissue groups
    tissue_df = aggregate_to_tissues(prop_df, marker_genes_convert)

    # Use timepoints only as a groupby key — never add to tissue_df as a column
    tp_fmt = [f"D{int(t):02d}" for t in time_points]
    tissue_df["tp_fmt"] = tp_fmt
    grouped = tissue_df.groupby("tp_fmt").mean()   # columns are tissue groups only

    target_cols = grouped.columns[grouped.mean() > 0.0]

    # 3. Plot-data CSV (tissue groups × timepoints — matches notebook output)
    plot_csv = os.path.join(OUT_DIR, f"Zhu_{tag}_plot_data.csv")
    grouped[target_cols].T.to_csv(plot_csv)
    print(f"    Saved: {plot_csv}")

    # 4. Stacked bar plot
    color_map = _build_color_map(target_cols)
    fig, ax = plt.subplots(figsize=(8, 4))
    grouped[target_cols].plot(kind="bar", stacked=True, ax=ax, color=color_map)
    ax.legend(loc="upper left", bbox_to_anchor=(1, 1), fontsize=7)
    ax.set_ylim(0, 1)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlabel("Day")
    ax.set_ylabel("Proportion")
    ax.set_title(f"Zhu et al. proteomics — {tag}")
    ax.grid(axis="y", alpha=0.4)
    plt.tight_layout()

    svg_path = os.path.join(OUT_DIR, f"Zhu_{tag}_proportion.svg")
    plt.savefig(svg_path, bbox_inches="tight")
    plt.close(fig)
    print(f"    Saved: {svg_path}")


# ── parameter grid ───────────────────────────────────────────────────────────
# Each entry: (algorithm_label, run_fn, kwargs_for_run_fn)
# run_fn signature: run_fn(bulk_matrix, reference_df, **kwargs)
#                   → list of pd.Series (one per sample)

# Minimum marker_score a gene must exceed to be included as a reference marker.
# The notebook uses 0.08 as the production value for proteomics.
MIN_MARKERSCORE_GRID = [0.0, 0.04, 0.08]

# Minimum number of expressed genes a cell type must have in the query data
# to be kept in the reference matrix.  Mirrors the cutoff explored in the
# notebook ("for proteomics we use min_marker=15").
MIN_MARKER_GRID = [0, 5, 10]

PARAM_GRID = [
    # ------------------------------------------------------------------ LSVR
    {
        "tag"   : "lsvr_C0.01_e0.1",
        "algo"  : "LSVR",
        "fn"    : run_multiprocess_deconvolution,
        "kwargs": {"param_C": 0.01,  "param_e": 0.1,  "max_workers": 4},
    },
    {
        "tag"   : "lsvr_C0.001_e0.1",
        "algo"  : "LSVR",
        "fn"    : run_multiprocess_deconvolution,
        "kwargs": {"param_C": 0.001, "param_e": 0.1,  "max_workers": 4},
    },
]


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    warnings.filterwarnings("ignore")

    # ── 1. Load Zhu et al. data (same as notebook) ──────────────────────────
    print("Loading Zhu et al. proteomics data …")
    os.chdir(os.path.join(ROOT, "notebooks"))   # utils.py uses relative paths
    target_df, time_points = zhu_ref_data()
    target_df = target_df.loc[:, ~target_df.isnull().values.all(axis=0)]
    target_df = target_df.astype(np.float32)
    os.chdir(ROOT)

    print(f"  Raw data: {target_df.shape[0]} samples × {target_df.shape[1]} genes")
    print(f"  Timepoints: {sorted(set(time_points))}")

    # ── 2. Load tissue-group mapping ────────────────────────────────────────
    marker_genes_convert = load_cell_to_group(
        filepath=os.path.join(DATA_DIR, "abbas_cell_markers.csv")
    )

    from scipy import stats as scipy_stats
    from itertools import permutations as perms

    # ── 3–4. Sweep both filter parameters, then algorithm configs ────────────
    # Outermost: min_markerscore — affects which genes qualify as markers
    # Middle:    min_marker      — affects which cell types survive (≥N genes)
    # Inner:     PARAM_GRID      — deconvolution algorithm / hyper-parameters
    summary_rows = []

    for min_markerscore in MIN_MARKERSCORE_GRID:
        for min_marker in MIN_MARKER_GRID:
            print(f"\n{'#'*60}")
            print(f"  min_markerscore = {min_markerscore}  |  min_marker = {min_marker}")

            # Both filter params change what goes into ref_mat / bulk_mat
            bulk_mat, ref_mat = build_aligned_matrices(
                target_df,
                min_markerscore=min_markerscore,
                min_marker=min_marker,
            )

            # Drop any non-gene columns (e.g. 'timepoints') before deconvolution
            non_gene_cols = [c for c in bulk_mat.columns if not str(c).startswith("WBGene")]
            if non_gene_cols:
                print(f"  Dropping non-gene columns from bulk matrix: {non_gene_cols}")
                bulk_mat = bulk_mat.drop(columns=non_gene_cols)

            n_marker_genes    = bulk_mat.shape[1]
            n_ref_cell_types  = ref_mat.shape[1]

            # Skip degenerate cases where no cell types survive the cutoffs
            if n_ref_cell_types == 0:
                print("  ⚠  No cell types remain — skipping.")
                continue

            for cfg in PARAM_GRID:
                algo = cfg["algo"]
                fn   = cfg["fn"]
                kw   = cfg["kwargs"]
                # Tag encodes markerscore, min_marker, and algorithm params
                tag  = f"ms{min_markerscore:.2f}_mm{min_marker:02d}_{cfg['tag']}"

                print(f"\n{'='*60}")
                print(f"  Running {algo} — {tag}")
                print(f"  Parameters: {kw}")

                results = fn(bulk_mat, ref_mat, **kw)

                # Assemble proportion DataFrame — cell-type columns only
                prop_df = pd.concat(results, axis=1).T
                prop_df.index = bulk_mat.index

                save_results(prop_df, time_points, tag, marker_genes_convert)

                # ── Summary statistics ───────────────────────────────────────
                tissue_df = aggregate_to_tissues(prop_df, marker_genes_convert)

                mean_by_tissue = tissue_df.mean()
                top3 = mean_by_tissue.nlargest(3)

                day1_mask = [t == min(time_points) for t in time_points]
                day1_data = tissue_df[day1_mask].values
                day1_corrs = []
                for i, j in list(perms(range(len(day1_data)), 2))[:50]:
                    r, _ = scipy_stats.spearmanr(day1_data[i], day1_data[j])
                    day1_corrs.append(r)
                mean_corr_d1 = np.mean(day1_corrs) if day1_corrs else float("nan")

                # Guard: fewer than 3 tissue groups when thresholds are strict
                def _t(i): return top3.index[i]  if len(top3) > i else "N/A"
                def _p(i): return round(float(top3.iloc[i]), 4) if len(top3) > i else 0.0

                summary_rows.append({
                    "min_markerscore" : min_markerscore,
                    "min_marker"      : min_marker,
                    "n_marker_genes"  : n_marker_genes,
                    "n_ref_cell_types": n_ref_cell_types,
                    "algorithm"       : algo,
                    "tag"             : tag,
                    "n_cell_types"    : len(prop_df.columns),
                    "top1_tissue"     : _t(0),
                    "top1_prop"       : _p(0),
                    "top2_tissue"     : _t(1),
                    "top2_prop"       : _p(1),
                    "top3_tissue"     : _t(2),
                    "top3_prop"       : _p(2),
                    "mean_within_day_corr_D1": round(mean_corr_d1, 4),
                })

                print(f"  Top tissues: {list(top3.index)}")
                print(f"  Day-1 within-day corr: {mean_corr_d1:.4f}")

    # ── 5. Save comparison summary ───────────────────────────────────────────
    summary_df = pd.DataFrame(summary_rows)
    summary_csv = os.path.join(OUT_DIR, "Zhu_algorithm_comparison.csv")
    summary_df.to_csv(summary_csv, index=False)
    print(f"\n{'='*60}")
    print(f"Summary saved: {summary_csv}")
    print()
    print(summary_df.to_string(index=False))

    # ── 6. Comparison plots — one SVG per sweep axis ─────────────────────────
    # Each SVG fixes the *other* parameter at its notebook default so that
    # the effect of each filter is visible independently.
    #   • Zhu_comparison_by_marker.svg     — fixes min_markerscore=0.08
    #   • Zhu_comparison_by_markerscore.svg — fixes min_marker=15
    # Kept as the canonical Zhu_algorithm_comparison.svg for backwards compat.

    from matplotlib.patches import Patch

    colors_algo = {"LSVR": "#4878CF",}

    def _facet_bar_plot(df, facet_col, facet_vals, strip_prefix_fmt, title, svg_path):
        """Draw one bar-chart panel per value in facet_vals."""
        n = len(facet_vals)
        fig, axes = plt.subplots(1, n, figsize=(4 * n, 4), sharey=True)
        if n == 1:
            axes = [axes]
        for ax, fval in zip(axes, facet_vals):
            sub = df[df[facet_col] == fval]
            if sub.empty:
                ax.set_visible(False)
                continue
            x_pos = range(len(sub))
            ax.bar(
                x_pos,
                sub["mean_within_day_corr_D1"],
                color=[colors_algo[a] for a in sub["algorithm"]],
                edgecolor="white", linewidth=0.8,
            )
            prefix = strip_prefix_fmt(fval)
            ax.set_xticks(list(x_pos))
            ax.set_xticklabels(
                sub["tag"].str.replace(prefix, "", regex=False),
                rotation=45, ha="right", fontsize=6,
            )
            # Panel title: show the facet value + resulting matrix dimensions
            n_genes = sub["n_marker_genes"].iloc[0]
            n_types = sub["n_ref_cell_types"].iloc[0]
            ax.set_title(
                f"{facet_col}={fval}\n({n_genes} genes, {n_types} cell types)",
                fontsize=8,
            )
            ax.set_ylim(0, 1)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.grid(axis="y", alpha=0.4)
        axes[0].set_ylabel("Mean within-day Spearman corr. (day 1)")
        fig.suptitle(title, fontsize=10)
        legend_elements = [Patch(facecolor=c, label=a) for a, c in colors_algo.items()]
        axes[-1].legend(handles=legend_elements, loc="lower right", fontsize=8)
        plt.tight_layout()
        plt.savefig(svg_path, bbox_inches="tight")
        plt.close(fig)
        print(f"Comparison plot: {svg_path}")

    # — Plot A: vary min_marker, fix min_markerscore at notebook default (0.08)
    df_ms_fixed = summary_df[summary_df["min_markerscore"] == 0.08]
    _facet_bar_plot(
        df_ms_fixed,
        facet_col="min_marker",
        facet_vals=MIN_MARKER_GRID,
        strip_prefix_fmt=lambda mm: f"ms0.08_mm{mm:02d}_",
        title="Zhu et al. proteomics — within-day agreement by min_marker (min_markerscore=0.08)",
        svg_path=os.path.join(OUT_DIR, "Zhu_comparison_by_marker.svg"),
    )

    # — Plot B: vary min_markerscore, fix min_marker at notebook default (15)
    df_mm_fixed = summary_df[summary_df["min_marker"] == 15]
    _facet_bar_plot(
        df_mm_fixed,
        facet_col="min_markerscore",
        facet_vals=MIN_MARKERSCORE_GRID,
        strip_prefix_fmt=lambda ms: f"ms{ms:.2f}_mm15_",
        title="Zhu et al. proteomics — within-day agreement by min_markerscore (min_marker=15)",
        svg_path=os.path.join(OUT_DIR, "Zhu_comparison_by_markerscore.svg"),
    )

    # Keep canonical name pointing at the min_marker plot (backwards compat)
    import shutil
    shutil.copy(
        os.path.join(OUT_DIR, "Zhu_comparison_by_marker.svg"),
        os.path.join(OUT_DIR, "Zhu_algorithm_comparison.svg"),
    )
    print(f"Canonical copy: {os.path.join(OUT_DIR, 'Zhu_algorithm_comparison.svg')}")

    print("\nAll done.")


if __name__ == "__main__":
    main()
