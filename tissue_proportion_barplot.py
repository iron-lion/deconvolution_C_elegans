"""
Bar plots of individual tissue proportions comparing:
  1. nuclei count fraction  (ground truth)
  2. cell volume fraction   (ground truth)
  3. swRNAseq N2 (ms=0.1, mm=0)    — best available; no per-algorithm breakdown
  4. Zhu      LSVR C=0.01 (ms=0.04, mm=5)  — best Spearman ρ_vol
"""
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT_DIR  = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(ROOT_DIR, 'results', 'Tissue_proportions')
DATA_DIR     = os.path.join(ROOT_DIR, 'data')
ZHU_DIR      = os.path.join(ROOT_DIR, 'results', 'Zhu')
EDER_DIR      = os.path.join(ROOT_DIR, 'results', 'Eder')
os.makedirs(OUT_DIR, exist_ok=True)

# ── Ground truth ───────────────────────────────────────────────────────────────
gt = pd.read_csv(os.path.join(DATA_DIR, 'Froehlich_tissue_volumnes.csv'), index_col=0)
gt.index = gt.index.str.strip()

GT_GROUPS = [
    'Intestine', 'Hypodermal cells', 'Muscle', 'Gonadal',
    'Uterine-vulval cells', 'Neurons', 'Spermatheca',
    'Spermatheca-Uterine junction', 'Vulval cell',
    'Intestinal-rectal valve', 'Glia', 'Distal tip cell', 'Coelomocytes',
]
vol_frac = gt.loc[GT_GROUPS, 'vol_fraction']
nuc_frac = gt.loc[GT_GROUPS, 'nuclei_fraction']

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

def aggregate_day1(path):
    df = pd.read_csv(path, index_col=0)
    d1 = df.iloc[:, 0]
    accum = {g: 0.0 for g in GT_GROUPS}
    for tissue, val in d1.items():
        cg = TISSUE_MAP.get(tissue.strip())
        if cg and cg in accum:
            accum[cg] += val
    return pd.Series(accum)

# ── Load the five series ───────────────────────────────────────────────────────
swrna = aggregate_day1(
    os.path.join(EDER_DIR, 'Eder_ms0.08_mm00_lsvr_C0.001_e0.1_plot_data.csv'))
zhu   = aggregate_day1(
    os.path.join(ZHU_DIR,     'Zhu_ms0.08_mm00_lsvr_C0.001_e0.1_plot_data.csv'))

SERIES = {
    'Nuclei fraction': nuc_frac,
    'Vol fraction':    vol_frac,
    'swRNAseq\n(ms=0.08, mm=0)':  swrna,
    'swProteomics\n(ms=0.08, mm=0)':     zhu,
}

COLORS = ['#1a9641', '#d7191c', '#ff7f00', '#1f78b4']

# ── Short tissue labels ────────────────────────────────────────────────────────
SHORT = {
    'Intestine':                    'Intestine',
    'Hypodermal cells':             'Hypodermis',
    'Muscle':                       'Muscle',
    'Gonadal':                      'Gonadal',
    'Uterine-vulval cells':         'Uterine-vulval',
    'Neurons':                      'Neurons',
    'Spermatheca':                  'Spermatheca',
    'Spermatheca-Uterine junction': 'Sp-Ut junction',
    'Vulval cell':                  'Vulval',
    'Intestinal-rectal valve':      'Int-rect valve',
    'Glia':                         'Glia',
    'Distal tip cell':              'Distal tip',
    'Coelomocytes':                 'Coelomocytes',
}
labels = [SHORT[g] for g in GT_GROUPS]

# ── Figure 1: grouped bar chart — all tissues together ────────────────────────
n_groups  = len(GT_GROUPS)
n_series  = len(SERIES)
bar_w     = 0.14
group_gap = 0.1
x = np.arange(n_groups) * (n_series * bar_w + group_gap)

fig, ax = plt.subplots(figsize=(22, 6))
for i, (name, series) in enumerate(SERIES.items()):
    vals = [series[g] for g in GT_GROUPS]
    offset = (i - n_series / 2 + 0.5) * bar_w
    bars = ax.bar(x + offset, vals, bar_w, label=name,
                  color=COLORS[i], alpha=0.88, edgecolor='white', linewidth=0.4)

ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=9, rotation=30, ha='right')
ax.set_ylabel('Proportion / Fraction', fontsize=11)
ax.set_title('Individual tissue proportions: ground truth vs. LSVR deconvolution (Day 1)',
             fontsize=13, fontweight='bold')
ax.legend(fontsize=8.5, ncol=n_series, loc='upper right',
          framealpha=0.9, edgecolor='gray')
ax.grid(axis='y', alpha=0.3, linestyle='--')
ax.set_ylim(0, max(max(s.values) for s in SERIES.values()) * 1.18)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
out1 = os.path.join(OUT_DIR, 'tissue_proportion_barplot.svg')
fig.savefig(out1, bbox_inches='tight')
plt.close()
print(f"Saved {out1}")

# ── Figure 2: one subplot per tissue ─────────────────────────────────────────
ncols = 5
nrows = -(-n_groups // ncols)   # ceiling division
fig, axes = plt.subplots(nrows, ncols,
                         figsize=(ncols * 3.2, nrows * 3.5),
                         sharey=False)
axes = axes.flatten()
fig.suptitle('Tissue-level proportion comparison: GT vs. LSVR deconvolution (Day 1)',
             fontsize=13, fontweight='bold', y=1.01)

x_ticks = np.arange(n_series)
short_names = [n.replace('\n', ' ') for n in SERIES.keys()]
tick_labels = [n.split('\n')[0] for n in SERIES.keys()]   # first line only

for idx, group in enumerate(GT_GROUPS):
    ax = axes[idx]
    vals = [SERIES[name][group] for name in SERIES]
    bars = ax.bar(x_ticks, vals, color=COLORS, alpha=0.88,
                  edgecolor='white', linewidth=0.5)
    for bar, v in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width() / 2, v + max(vals) * 0.02,
                f'{v:.3f}', ha='center', va='bottom', fontsize=6.5, rotation=0)
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(tick_labels, fontsize=7, rotation=30, ha='right')
    ax.set_title(SHORT[group], fontsize=10, fontweight='bold')
    ax.set_ylabel('Proportion', fontsize=8)
    ax.set_ylim(0, max(vals) * 1.28 if max(vals) > 0 else 0.01)
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

# hide empty subplots
for idx in range(n_groups, len(axes)):
    axes[idx].set_visible(False)

# shared legend below figure
legend_handles = [
    plt.Rectangle((0, 0), 1, 1, color=COLORS[i], alpha=0.88)
    for i in range(n_series)
]
fig.legend(legend_handles, short_names,
           loc='lower center', ncol=n_series, fontsize=8.5,
           bbox_to_anchor=(0.5, -0.04), framealpha=0.9, edgecolor='gray')

plt.tight_layout()
out2 = os.path.join(OUT_DIR, 'tissue_proportion_per_tissue.svg')
fig.savefig(out2, bbox_inches='tight')
plt.close()
print(f"Saved {out2}")

# ── Table: all values ─────────────────────────────────────────────────────────
rows = []
for g in GT_GROUPS:
    row = {'tissue': SHORT[g]}
    for name, series in SERIES.items():
        col = name.replace('\n', ' ')
        row[col] = round(series[g], 4)
    rows.append(row)
df_table = pd.DataFrame(rows)
out3 = os.path.join(OUT_DIR, 'tissue_proportion_table.csv')
df_table.to_csv(out3, index=False)
print(f"Saved {out3}")
print(df_table.to_string(index=False))
