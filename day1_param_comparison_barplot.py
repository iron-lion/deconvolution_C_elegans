"""
day1_param_comparison_barplot.py
=================================
Day-1 raw LSVR Spearman ρ, Pearson r, and Jensen-Shannon divergence vs.
two ground truths for all available parameter cutoff combinations across
Zhu, Eder.

Ground truths (both from Froehlich_tissue_volumnes.csv — consistent with
combined_ground_truth_comparison.py):
  vol_fraction     — cell volume fraction 
  nuclei_fraction  — nuclei count fraction

Parameter cutoffs shown on x-axis:
  ms  = marker score threshold (0.00, 0.04, 0.08)
  mm  = minimum markers per cell type (0, 5, 10)
  C   = SVR regularisation (0.001, 0.01)

JS divergence notes:
  - Computed as JS distance (sqrt of JSD), bounded [0, 1]
  - Both prediction and GT vectors normalised to sum to 1 before computation
  - Lower JS distance = more similar distributions (opposite of correlation)

Outputs to results/cross_dataset/:
  day1_param_all.svg        ← 3×2 combined (Spearman+Pearson+JS × vol+nuc)
  day1_param_rho_vol.svg    ← Spearman ρ vs. cell volume
  day1_param_r_vol.svg      ← Pearson r vs. cell volume
  day1_param_rho_nuc.svg    ← Spearman ρ vs. nuclei count
  day1_param_r_nuc.svg      ← Pearson r vs. nuclei count
  day1_param_js_vol.svg     ← JS distance vs. cell volume
  day1_param_js_nuc.svg     ← JS distance vs. nuclei count
  day1_param_metrics.csv    ← full table
"""

import os, sys, re
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats
from scipy.spatial.distance import jensenshannon

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(ROOT_DIR, "src"))
DATA_DIR    = os.path.join(ROOT_DIR, "data")
OUT_DIR     = os.path.join(ROOT_DIR, "results", "ground_truth")
os.makedirs(OUT_DIR, exist_ok=True)

# ── Ground truth: same source as combined_ground_truth_comparison.py ──────────
GT_GROUPS = [
    'Intestine', 'Hypodermal cells', 'Muscle', 'Gonadal',
    'Uterine-vulval cells', 'Neurons', 'Spermatheca',
    'Spermatheca-Uterine junction', 'Vulval cell',
    'Intestinal-rectal valve', 'Glia', 'Distal tip cell', 'Coelomocytes',
]

gt = pd.read_csv(os.path.join(DATA_DIR, "Froehlich_tissue_volumnes.csv"), index_col=0)
gt.index = gt.index.str.strip()
vol_frac = gt.loc[GT_GROUPS, "vol_fraction"]   # raw, sum ≈ 0.98
nuc_frac = gt.loc[GT_GROUPS, "nuclei_fraction"] # raw, from cell counts

# ── Tissue aggregation (same TISSUE_MAP as combined_ground_truth_comparison.py)
TISSUE_MAP = {
    'Neurons':                      'Neurons',
    'Muscle':                       'Muscle',
    'Anal muscle':                  'Muscle',
    'Hypodermal cells':             'Hypodermal cells',
    'Intestine':                    'Intestine',
    'Intestinal-rectal valve':      'Intestinal-rectal valve',
    'Gonadal':                      'Gonadal',
    'Distal tip cell':              'Distal tip cell',
    'Spermatheca':                  'Spermatheca',
    'Spermatheca-Uterine junction': 'Spermatheca-Uterine junction',
    'Uterine-vulval cells':         'Uterine-vulval cells',
    'Vulval cell':                  'Vulval cell',
    'Glia':                         'Glia',
    'Amphid socket':                'Glia',
    'Phasmid socket':               'Glia',
    'Unassigned sheath cells':      'Glia',
    'cephalic sheath  and unidentified glial': 'Glia',
    'Coelomocytes':                 'Coelomocytes',
    'Excretory cell':               'Coelomocytes',
}

def aggregate_to_groundtruth(series):
    """Collapse plot_data cell-group rows into 13 GT groups."""
    accum = {g: 0.0 for g in GT_GROUPS}
    for tissue, val in series.items():
        tg = TISSUE_MAP.get(str(tissue).strip())
        if tg and tg in accum:
            accum[tg] += float(val)
    return pd.Series(accum, index=GT_GROUPS)

# ── File discovery ────────────────────────────────────────────────────────────
_PAT = re.compile(r"ms([\d.]+)_mm(\d+)_lsvr_C([\d.]+)_e")

SEARCH_DIRS = {
    "Zhu":     (os.path.join(ROOT_DIR, "results", "Zhu"),          0),
    "Eder": (os.path.join(ROOT_DIR, "results", "Eder"), 0),
}

# ── Compute correlations + JS divergence ─────────────────────────────────────
v = vol_frac.values
n = nuc_frac.values

def js_dist(p, q, eps=1e-10):
    """JS distance (sqrt JSD) after normalising both vectors to sum to 1."""
    p = np.array(p, dtype=float) + eps
    q = np.array(q, dtype=float) + eps
    return float(jensenshannon(p / p.sum(), q / q.sum()))

rows = []

for ds, (folder, _) in SEARCH_DIRS.items():
    for fname in sorted(os.listdir(folder)):
        if not fname.endswith("_plot_data.csv"):
            continue
        m = _PAT.search(fname)
        if not m:
            continue
        ms, mm, C = m.group(1), m.group(2), m.group(3)
        param_tag = f"ms{ms}_mm{mm}_C{C}"

        df_file = pd.read_csv(os.path.join(folder, fname), index_col=0)
        d1  = df_file.iloc[:, 0]
        agg = aggregate_to_groundtruth(d1)
        p   = agg.values

        rows.append({
            "dataset":  ds,
            "param_tag": param_tag,
            "rho_vol":  stats.spearmanr(p, v)[0],
            "r_vol":    stats.pearsonr(p,  v)[0],
            "rho_nuc":  stats.spearmanr(p, n)[0],
            "r_nuc":    stats.pearsonr(p,  n)[0],
            "js_vol":   js_dist(p, v),
            "js_nuc":   js_dist(p, n),
        })

df = pd.DataFrame(rows)
df.to_csv(os.path.join(OUT_DIR, "day1_param_metrics.csv"), index=False)
print(f"Collected {len(df)} rows across {df['dataset'].nunique()} datasets.")
print(df[["dataset", "param_tag", "rho_vol", "r_vol", "rho_nuc", "r_nuc",
          "js_vol", "js_nuc"]].to_string(index=False))

# ── Layout constants ──────────────────────────────────────────────────────────
DATASET_ORDER  = ["Eder", "Zhu"]
DATASET_LABELS = {"Zhu": "Zhu (LFQ)",
                  "Eder": "Eder (swRNAseq)"}
present_ds = [d for d in DATASET_ORDER if d in df["dataset"].values]
all_tags   = sorted(df["param_tag"].unique())

DS_COLORS = {
    "Zhu":      "#1f78b4",
    "Eder": "#ff7f00",
}

# ── Plotting helper ───────────────────────────────────────────────────────────
def param_barplot(ax, metric, title, ylabel, ylim=(-0.2, 1.05)):
    tags_present = sorted(df["param_tag"].unique())
    n_tags = len(tags_present)
    n_ds   = len(present_ds)
    width  = 0.8 / n_ds
    x      = np.arange(n_tags)

    for i, ds in enumerate(present_ds):
        sub  = df[df["dataset"] == ds].set_index("param_tag")
        vals = [sub.loc[t, metric] if t in sub.index else np.nan
                for t in tags_present]
        offset = (i - n_ds / 2 + 0.5) * width
        bars = ax.bar(x + offset, vals, width * 0.88,
                      color=DS_COLORS[ds], edgecolor="black",
                      linewidth=0.5, alpha=0.88, label=DATASET_LABELS[ds])
        for bar, val in zip(bars, vals):
            if np.isfinite(val):
                ax.text(bar.get_x() + bar.get_width() / 2,
                        val + 0.012, f"{val:.2f}",
                        ha="center", va="bottom", fontsize=5.5,
                        color="black", fontweight="bold", rotation=90)

    ax.set_xticks(x)
    ax.set_xticklabels(tags_present, fontsize=7.5, rotation=40, ha="right")
    ax.set_ylabel(ylabel, fontsize=10)
    ax.set_title(title, fontsize=10.5, fontweight="bold", pad=5)
    ax.set_ylim(*ylim)
    ax.axhline(0, color="black", lw=0.6, alpha=0.4)
    ax.grid(axis="y", alpha=0.25, linestyle="--")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

legend_handles = [mpatches.Patch(facecolor=DS_COLORS[d], edgecolor="black",
                                  label=DATASET_LABELS[d]) for d in present_ds]

# JS distance: lower = better → invert y-axis or annotate clearly
def param_barplot_js(ax, metric, title, ylabel, ylim=(0, 0.55)):
    """Like param_barplot but annotated for JS distance (lower = better)."""
    tags_present = sorted(df["param_tag"].unique())
    n_tags = len(tags_present)
    n_ds   = len(present_ds)
    width  = 0.8 / n_ds
    x      = np.arange(n_tags)

    for i, ds in enumerate(present_ds):
        sub  = df[df["dataset"] == ds].set_index("param_tag")
        vals = [sub.loc[t, metric] if t in sub.index else np.nan
                for t in tags_present]
        offset = (i - n_ds / 2 + 0.5) * width
        bars = ax.bar(x + offset, vals, width * 0.88,
                      color=DS_COLORS[ds], edgecolor="black",
                      linewidth=0.5, alpha=0.88, label=DATASET_LABELS[ds])
        for bar, val in zip(bars, vals):
            if np.isfinite(val):
                ax.text(bar.get_x() + bar.get_width() / 2,
                        val + 0.005, f"{val:.2f}",
                        ha="center", va="bottom", fontsize=5.5,
                        color="black", fontweight="bold", rotation=90)

    ax.set_xticks(x)
    ax.set_xticklabels(tags_present, fontsize=7.5, rotation=40, ha="right")
    ax.set_ylabel(ylabel + "  (↓ better)", fontsize=10)
    ax.set_title(title, fontsize=10.5, fontweight="bold", pad=5)
    ax.set_ylim(*ylim)
    ax.grid(axis="y", alpha=0.25, linestyle="--")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

# ── 3×2 combined figure ───────────────────────────────────────────────────────
fig, axes = plt.subplots(3, 2, figsize=(22, 19))
fig.suptitle(
    "Day-1 Raw LSVR — Correlation & JS Divergence across Parameter Cutoffs\n"
    "GT: cell volume fraction & nuclei count fraction (tissue volume data)\n"
    "ms = marker score  |  mm = min markers per cell type  |  C = SVR regularisation\n"
    "JS distance: lower = better  |  Spearman/Pearson: higher = better",
    fontsize=12, fontweight="bold",
)

param_barplot(axes[0, 0], "rho_vol",
              "Spearman ρ  vs.  Cell Volume Fraction",
              "Spearman ρ")
param_barplot(axes[0, 1], "rho_nuc",
              "Spearman ρ  vs.  Nuclei Count Fraction",
              "Spearman ρ")
param_barplot(axes[1, 0], "r_vol",
              "Pearson r  vs.  Cell Volume Fraction",
              "Pearson r")
param_barplot(axes[1, 1], "r_nuc",
              "Pearson r  vs.  Nuclei Count Fraction",
              "Pearson r")
param_barplot_js(axes[2, 0], "js_vol",
                 "JS Distance  vs.  Cell Volume Fraction",
                 "JS Distance")
param_barplot_js(axes[2, 1], "js_nuc",
                 "JS Distance  vs.  Nuclei Count Fraction",
                 "JS Distance")

fig.legend(handles=legend_handles, loc="lower center",
           ncol=len(present_ds), fontsize=10,
           bbox_to_anchor=(0.5, -0.01), framealpha=0.92, edgecolor="gray")
plt.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "day1_param_all.svg"), bbox_inches="tight")
plt.close()
print("\nSaved day1_param_all.svg")

# ── Individual figures ────────────────────────────────────────────────────────
corr_panels = [
    ("rho_vol", "Spearman ρ vs. Cell Volume Fraction",  "day1_param_rho_vol.svg", False),
    ("r_vol",   "Pearson r vs. Cell Volume Fraction",   "day1_param_r_vol.svg",   False),
    ("rho_nuc", "Spearman ρ vs. Nuclei Count Fraction", "day1_param_rho_nuc.svg", False),
    ("r_nuc",   "Pearson r vs. Nuclei Count Fraction",  "day1_param_r_nuc.svg",   False),
    ("js_vol",  "JS Distance vs. Cell Volume Fraction", "day1_param_js_vol.svg",  True),
    ("js_nuc",  "JS Distance vs. Nuclei Count Fraction","day1_param_js_nuc.svg",  True),
]

for metric, label, fname, is_js in corr_panels:
    fig, ax = plt.subplots(figsize=(max(14, len(all_tags) * 1.3), 6.5))
    fig.suptitle(
        f"Day-1 Raw LSVR — {label}\n"
        "Across All Parameter Cutoffs and Datasets\n"
        "GT from tissue volume data ",
        fontsize=11, fontweight="bold",
    )
    if is_js:
        param_barplot_js(ax, metric, "", label.split(" vs.")[0])
    else:
        param_barplot(ax, metric, "", label.split(" vs.")[0])
    ax.set_xlabel("Parameter cutoff  (ms = marker score, mm = min markers, C = regularisation)",
                  fontsize=9)
    fig.legend(handles=legend_handles, loc="lower center",
               ncol=len(present_ds), fontsize=9,
               bbox_to_anchor=(0.5, -0.06), framealpha=0.92, edgecolor="gray")
    plt.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, fname), bbox_inches="tight")
    plt.close()
    print(f"Saved {fname}")

print(f"\nAll done. Outputs in {OUT_DIR}")
