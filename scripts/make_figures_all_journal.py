#!/usr/bin/env python3
"""Publication figures for the all-241 prophage analysis (journal style)."""
import json
import os
import re
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import Phylo
from matplotlib import font_manager as fm
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from scipy import stats
from scipy.spatial.distance import pdist, squareform

BASE = "/home/bioinfo/cdprophage"
RESULT = os.path.join(BASE, "scaleup_241", "results")
OUT = os.path.join(BASE, "scaleup_241", "figures")
os.makedirs(OUT, exist_ok=True)

# Arial is the journal-preferred sans-serif face; register the Windows copy.
for _f in ("arial.ttf", "arialbd.ttf", "ariali.ttf", "arialbi.ttf"):
    _p = os.path.join("/mnt/c/Windows/Fonts", _f)
    if os.path.exists(_p):
        fm.fontManager.addfont(_p)

plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 8,
    "axes.labelsize": 8,
    "axes.titlesize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 150,
    "savefig.dpi": 600,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "mathtext.fontset": "custom",
    "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold",
    "mathtext.default": "it",
})

CLADE_COLORS = {"1": "#1f77b4", "2": "#ff7f0e", "3": "#2ca02c", "4": "#d62728", "5": "#8c564b"}
CONF_COLORS = {"high": "#1f4e79", "medium": "#5b9bd5", "low": "#c9c9c9"}
PRESENT_COLOR = "#08306b"
W = 6.22  # final print width (15.8 cm) in inches


def tag(ax, letter):
    ax.set_title("(%s)" % letter, loc="left", fontsize=10, fontweight="bold", pad=4)


def fmt_p(p):
    """Format a P value for a panel annotation without exponent notation."""
    if p >= 0.001:
        return "$P$ = %s" % ("%.3g" % p)
    exp = int(np.floor(np.log10(p)))
    mant = p / 10 ** exp
    if mant >= 9.95:
        mant, exp = mant / 10, exp + 1
    return "$P$ = %.1f $\\times$ 10$^{%d}$" % (mant, exp)


def save(fig, name):
    fig.savefig(os.path.join(OUT, name + ".png"), dpi=600, bbox_inches="tight", pad_inches=0.02)
    fig.savefig(os.path.join(OUT, name + ".pdf"), bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)


genomes = pd.read_csv(os.path.join(RESULT, "genome_summary_analysis.csv"), dtype={"ST": str, "mlst_clade": str})
details = pd.read_csv(os.path.join(RESULT, "prophage_details.csv"), dtype={"ST": str, "mlst_clade": str})

# ---------- Fig 1 ----------
fig, axes = plt.subplots(1, 2, figsize=(W, 2.95), gridspec_kw={"width_ratios": [1.35, 1]})
ax = axes[0]
clades = sorted(genomes["mlst_clade"].unique(), key=int)
data = [genomes[genomes["mlst_clade"] == c]["n_prophages"].values for c in clades]
ax.boxplot(data, positions=range(len(clades)), widths=0.55, patch_artist=True,
           showfliers=False, medianprops={"color": "#222222", "linewidth": 1.2},
           whiskerprops={"color": "#666666"}, capprops={"color": "#666666"},
           boxprops={"facecolor": "#e8eef5", "edgecolor": "#666666", "linewidth": 0.8})
rng = np.random.default_rng(3)
for i, (c, vals) in enumerate(zip(clades, data)):
    x = rng.normal(i, 0.075, len(vals))
    ax.scatter(x, vals, s=7, alpha=0.65, color=CLADE_COLORS[c], edgecolor="white", linewidth=0.2, zorder=3)
ax.set_xticks(range(len(clades)))
ax.set_xticklabels([f"Clade {c}\n($n$={len(genomes[genomes['mlst_clade']==c])})" for c in clades])
ax.set_ylabel("Prophage regions per genome")
tag(ax, "a")

ax = axes[1]
metrics = [
    ("Any prophage", genomes["carriage"].mean()),
    ("High-confidence", genomes["high_carriage"].mean()),
    ("Complete-like", genomes["complete_carriage"].mean()),
    ("Reference-typeable", (genomes["n_typeable"] > 0).mean()),
]
labels = [m[0] for m in metrics]
vals = [m[1] * 100 for m in metrics]
ax.barh(range(len(vals))[::-1], vals, color=["#1f4e79", "#2e75b6", "#5b9bd5", "#a6a6a6"], height=0.6)
for y, v in zip(range(len(vals))[::-1], vals):
    ax.text(v + 1.5, y, f"{v:.1f}%", va="center", fontsize=8)
ax.set_yticks(range(len(vals))[::-1])
ax.set_yticklabels(labels)
ax.set_xlim(0, 112)
ax.set_xlabel("Genomes carrying \u22651 region (%)")
tag(ax, "b")
fig.tight_layout()
save(fig, "fig1_carriage_by_clade")

# ---------- Fig 2 ----------
fig, axes = plt.subplots(1, 2, figsize=(W, 2.8))
ax = axes[0]
bins = np.arange(0, max(80000, details["length_bp"].max()) + 5000, 5000)
bottom = np.zeros(len(bins) - 1)
for conf in ["high", "medium", "low"]:
    vals = details[details["confidence"] == conf]["length_bp"].values
    hist, _ = np.histogram(vals, bins=bins)
    ax.bar(bins[:-1] + 2500, hist, width=4500, bottom=bottom, color=CONF_COLORS[conf],
           label=conf.capitalize(), edgecolor="white", linewidth=0.4)
    bottom += hist
ax.set_xlabel("Prophage region length (kb)")
ax.set_ylabel("Number of regions")
ax.set_xticks(np.arange(0, bins[-1] + 1, 10000))
ax.set_xticklabels([f"{int(x/1000)}" for x in np.arange(0, bins[-1] + 1, 10000)])
ax.legend(frameon=False, title="Confidence")
tag(ax, "a")

ax = axes[1]
for conf in ["high", "medium", "low"]:
    sub = details[details["confidence"] == conf]
    ax.scatter(sub["length_bp"] / 1000, sub["gc_percent"], s=8, color=CONF_COLORS[conf],
               edgecolor="white", linewidth=0.2, label=conf.capitalize(), alpha=0.75)
mean_gc = genomes["genome_gc_percent"].mean()
ax.axhline(mean_gc, color="#d62728", linestyle="--", linewidth=1)
ax.text(96, mean_gc + 0.7, f"Mean host GC ({mean_gc:.1f}%)", color="#d62728", fontsize=8, ha="right",
        bbox={"facecolor": "white", "edgecolor": "none", "pad": 1})
ax.set_xlabel("Prophage region length (kb)")
ax.set_ylabel("Prophage GC (%)")
tag(ax, "b")
fig.tight_layout()
save(fig, "fig2_size_gc")

# ---------- Fig 3 ----------
clusters = sorted(details["cluster_id"].dropna().unique())
presence = pd.DataFrame(0, index=genomes["assembly_accession"], columns=clusters)
for _, r in details.iterrows():
    presence.loc[r["assembly_accession"], r["cluster_id"]] = 1
counts = presence.sum(axis=0).sort_values(ascending=False)
sharing = Counter()
for v in counts:
    sharing["1" if v == 1 else ("2-3" if v <= 3 else ("4-10" if v <= 10 else ">10"))] += 1

fig, axes = plt.subplots(1, 3, figsize=(W, 2.45), gridspec_kw={"width_ratios": [1.1, 0.9, 1.2]})
ax = axes[0]
ax.bar(range(1, len(counts) + 1), counts.values, color="#2e75b6", width=0.9)
ax.set_xlabel("Cluster rank")
ax.set_ylabel("Genomes carrying cluster")
tag(ax, "a")

ax = axes[1]
order = ["1", "2-3", "4-10", ">10"]
vals = [sharing.get(k, 0) for k in order]
ax.bar(range(len(order)), vals, color=["#c9c9c9", "#9dc3e6", "#2e75b6", "#1f4e79"], width=0.65)
for i, v in enumerate(vals):
    ax.text(i, v + max(vals) * 0.03, str(v), ha="center", fontsize=8)
ax.set_xticks(range(len(order)))
ax.set_xticklabels(["1", "2-3", "4-10", ">10"])
ax.set_ylabel("Number of clusters")
ax.set_ylim(0, max(vals) * 1.18)
tag(ax, "b")

ax = axes[2]
acc = pd.read_csv(os.path.join(RESULT, "cluster_accumulation.csv"))
ax.plot(acc["n_genomes"], acc["clusters_mean"], color="#1f4e79", linewidth=1.8)
ax.fill_between(acc["n_genomes"], acc["clusters_p05"], acc["clusters_p95"], color="#9dc3e6", alpha=0.35)
ax.set_xlabel("Number of genomes sampled")
ax.set_ylabel("Cumulative clusters")
ax.axhline(len(counts), color="#c0392b", linestyle="--", linewidth=1)
ax.text(3, len(counts) + 3, f"{len(counts)} clusters observed", color="#c0392b", fontsize=8, va="bottom")
ax.set_ylim(0, len(counts) * 1.12)
tag(ax, "c")
fig.tight_layout()
save(fig, "fig3_cluster_sharing")

# ---------- Fig 4 ----------
tree = Phylo.read(os.path.join(RESULT, "mash_nj_tree.nwk"), "newick")
tree_order = [t.name for t in tree.get_terminals() if t.name in set(genomes["assembly_accession"])]
missing = [a for a in genomes["assembly_accession"] if a not in tree_order]
tree_order += missing
meta_all = genomes.set_index("assembly_accession")
# group genomes by clade, preserving Mash-tree order within each clade
order = []
for c in clades:
    order += [a for a in tree_order if str(meta_all.loc[a, "mlst_clade"]) == c]
top = list(counts.index[:60])
# rows = clusters, columns = genomes: this keeps every cluster label horizontal
mat = presence.loc[order, top].T.values
clade_codes = np.array([int(meta_all.loc[a, "mlst_clade"]) for a in order])
bounds = np.cumsum([sum(1 for a in order if str(meta_all.loc[a, "mlst_clade"]) == c) for c in clades])

fig = plt.figure(figsize=(W, 7.8))
gs = fig.add_gridspec(2, 2, width_ratios=[1, 0.16], height_ratios=[0.055, 1],
                      hspace=0.02, wspace=0.03)
ax_top = fig.add_subplot(gs[0, 0])
ax = fig.add_subplot(gs[1, 0], sharex=ax_top)
ax_right = fig.add_subplot(gs[1, 1])

ax_top.imshow(clade_codes.reshape(1, -1), aspect="auto", interpolation="nearest",
              cmap=ListedColormap([CLADE_COLORS[str(i)] for i in range(1, 6)]), vmin=1, vmax=5)
ax_top.set_yticks([])
ax_top.set_xticks([])
for spine in ax_top.spines.values():
    spine.set_visible(False)
for c, start, end in zip(clades, np.r_[0, bounds[:-1]], bounds):
    ax_top.text((start + end - 1) / 2, 0, str(c), ha="center", va="center",
                color="white", fontsize=7, fontweight="bold")

ax.imshow(mat, aspect="auto", interpolation="nearest",
          cmap=ListedColormap(["#f2f2f2", PRESENT_COLOR]), vmin=0, vmax=1)
ax.set_yticks(range(len(top)))
ax.set_yticklabels(top, fontsize=8)
ax.set_xticks([])
ax.set_xlabel("Genomes ($n$ = 241), grouped by MLST clade")
ax.set_ylabel("Prophage cluster (ordered by prevalence)")
for spine in ax.spines.values():
    spine.set_linewidth(0.4)
for b in bounds[:-1]:
    ax.axvline(b - 0.5, color="white", linewidth=0.9)
    ax_top.axvline(b - 0.5, color="white", linewidth=0.9)

prev = counts.loc[top].values
ax_right.barh(range(len(top)), prev, color="#4292c6", height=0.78)
ax_right.set_xlim(0, max(prev) * 1.32)
ax_right.set_ylim(len(top) - 0.5, -0.5)
ax_right.set_yticks([])
ax_right.set_xticks([0, 20, 40])
ax_right.set_xlabel("Genomes", fontsize=8)
ax_right.tick_params(axis="x", labelsize=7)
for y, v in enumerate(prev):
    ax_right.text(v + max(prev) * 0.03, y, str(int(v)), va="center", fontsize=7)
for s in ("top", "right"):
    ax_right.spines[s].set_visible(False)
ax_right.set_title("Prevalence", fontsize=8, loc="left", pad=2)

legend_elements = [Patch(facecolor=CLADE_COLORS[str(i)], label=f"Clade {i}") for i in range(1, 6)]
fig.legend(handles=legend_elements, loc="lower center", ncol=5, frameon=False,
           bbox_to_anchor=(0.45, -0.005), handlelength=1.0, columnspacing=1.2)
fig.subplots_adjust(bottom=0.07, right=0.97)
save(fig, "fig4_cluster_heatmap")

# ---------- Fig 5 ----------
gc = pd.read_csv(os.path.join(RESULT, "prophage_gene_content.csv"), dtype={"ST": str, "mlst_clade": str})
fam = {}
cur = None
for line in open(os.path.join(RESULT, "prophage_proteins_nr.fasta.clstr")):
    line = line.rstrip("\n")
    if line.startswith(">Cluster"):
        cur = int(line.split()[1])
        continue
    m = re.search(r">([^|]+)\|", line)
    if m:
        fam.setdefault(m.group(1), set()).add("GF%05d" % cur)
pids = sorted(fam)
allgf = sorted({g for s in fam.values() for g in s})
gindex = {g: i for i, g in enumerate(allgf)}
X = np.zeros((len(pids), len(allgf)), dtype=np.int8)
for i, p in enumerate(pids):
    for g in fam[p]:
        X[i, gindex[g]] = 1
gc = gc.set_index("prophage_id").loc[pids].reset_index()
D = squareform(pdist(X, metric="jaccard"))
n = D.shape[0]
A = -0.5 * D ** 2
Jc = np.eye(n) - np.ones((n, n)) / n
B = Jc @ A @ Jc
vals, vecs = np.linalg.eigh(B)
order = np.argsort(vals)[::-1]
vals = vals[order]
vecs = vecs[:, order]
coords = vecs[:, :2] * np.sqrt(np.maximum(vals[:2], 0))
explained = np.maximum(vals[:2], 0) / vals[vals > 0].sum()
json.dump({"pcoa1_pct": round(float(explained[0] * 100), 2),
           "pcoa2_pct": round(float(explained[1] * 100), 2),
           "regions": int(n)},
          open(os.path.join(RESULT, "pcoa_variance.json"), "w"), indent=2)

meta = gc.set_index("prophage_id").loc[gc["prophage_id"]]
fig = plt.figure(figsize=(W, 5.9))
gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1], hspace=0.45, wspace=0.3)
ax = fig.add_subplot(gs[0, 0])
for c in clades:
    sel = meta["mlst_clade"].astype(str) == c
    ax.scatter(coords[sel, 0], coords[sel, 1], s=6, color=CLADE_COLORS[c], alpha=0.7,
               edgecolor="white", linewidth=0.15, label=f"Clade {c} ($n$={sel.sum()})")
ax.set_xlabel(f"PCoA1 ({explained[0]*100:.1f}% of variance)")
ax.set_ylabel(f"PCoA2 ({explained[1]*100:.1f}% of variance)")
ax.text(0.02, 0.02, f"Axes 1-2: {100*(explained[0]+explained[1]):.1f}%",
        transform=ax.transAxes, ha="left", va="bottom", fontsize=8)
ax.legend(frameon=False, markerscale=1.5)
tag(ax, "a")

ax = fig.add_subplot(gs[0, 1])
clade_arr = meta["mlst_clade"].astype(str).to_numpy()
within, between = [], []
for i in range(n):
    for j in range(i + 1, n):
        sim = 1.0 - D[i, j]
        (within if clade_arr[i] == clade_arr[j] else between).append(sim)
uw, pw = stats.mannwhitneyu(within, between, alternative="greater")
ax.violinplot([within, between], positions=[0, 1], widths=0.7, showextrema=False)
ax.boxplot([within, between], positions=[0, 1], widths=0.18, patch_artist=True,
           showfliers=False, medianprops={"color": "#222222", "linewidth": 1},
           boxprops={"facecolor": "white", "edgecolor": "#444444"},
           whiskerprops={"color": "#444444"}, capprops={"color": "#444444"})
ax.set_xticks([0, 1])
ax.set_xticklabels(["Within\nclade", "Between\nclades"])
ax.set_ylabel("Pairwise Jaccard similarity")
ax.text(0.03, 0.95, f"mean {np.mean(within):.3f} vs {np.mean(between):.3f}\nMann-Whitney {fmt_p(pw)}",
        transform=ax.transAxes, fontsize=8, va="top")
json.dump({"within_mean": round(float(np.mean(within)), 4),
           "between_mean": round(float(np.mean(between)), 4),
           "n_pairs_within": int(len(within)), "n_pairs_between": int(len(between)),
           "mannwhitney_p": float(pw)},
          open(os.path.join(RESULT, "gene_sharing_stats.json"), "w"), indent=2)
tag(ax, "b")

ax = fig.add_subplot(gs[1, :])
pn = json.load(open(os.path.join(RESULT, "protein_network_summary.json")))
cats = ["phage_core", "defence", "MGE", "virulence", "AMR"]
labels = {"phage_core": "Phage core", "defence": "Defence", "MGE": "MGE",
          "virulence": "Virulence", "AMR": "AMR"}
tot = pn["gene_families_by_category"]
cross = pn["cross_clade_gene_families"]
vals2 = [100 * cross[c] / tot[c] if tot[c] else 0 for c in cats]
colors = ["#1f4e79", "#c0392b", "#7f8c8d", "#e67e22", "#2e86c1"]
ypos = np.arange(len(cats))[::-1]
ax.barh(ypos, vals2, color=colors, height=0.6)
for y, v, c in zip(ypos, vals2, cats):
    ax.text(v + 1.0, y, f"{v:.0f}% ({cross[c]}/{tot[c]})", va="center", fontsize=8)
ax.set_yticks(ypos)
ax.set_yticklabels([labels[c] for c in cats])
ax.set_xlabel("% of gene families in ≥2 clades")
ax.set_xlim(0, max(60, max(vals2) * 1.35))
tag(ax, "c")
save(fig, "fig5_protein_network")

# ---------- Fig 6 ----------
cargo = pd.read_csv(os.path.join(RESULT, "prophage_cargo_summary.csv"), dtype={"ST": str, "mlst_clade": str})
cs = json.load(open(os.path.join(RESULT, "cargo_summary.json")))
fig = plt.figure(figsize=(W, 6.0))
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.05], hspace=0.5, wspace=0.3)
ax = fig.add_subplot(gs[0, 0])
cats = ["defence", "virulence", "AMR", "MGE"]
labels = ["Defence", "Virulence", "AMR", "MGE"]
colmap = {"defence": "has_defence", "virulence": "has_virulence", "AMR": "has_amr", "MGE": "has_mge"}
high_vals = [100 * cargo[cargo["confidence"] == "high"][colmap[c]].mean() for c in cats]
all_vals = [100 * cargo[colmap[c]].mean() for c in cats]
x = np.arange(len(cats))
ax.bar(x - 0.19, all_vals, width=0.36, color="#9dc3e6", label=f"All regions ($n$={len(cargo)})")
ax.bar(x + 0.19, high_vals, width=0.36, color="#1f4e79", label=f"High-confidence ($n$={(cargo['confidence']=='high').sum()})")
for i, (a, h) in enumerate(zip(all_vals, high_vals)):
    ax.text(i - 0.19, a + 1.5, f"{a:.0f}", ha="center", fontsize=8)
    ax.text(i + 0.19, h + 1.5, f"{h:.0f}", ha="center", fontsize=8)
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.set_ylabel("% of prophage regions")
ax.set_ylim(0, 108)
ax.legend(frameon=False, loc="upper right", fontsize=8)
tag(ax, "a")

ax = fig.add_subplot(gs[0, 1])
en = cs["enrichment_vs_chromosome"]
en_map = {"defence": "defence", "virulence": "virulence", "AMR": "amr", "MGE": "mge"}
x = np.arange(len(cats))
pp = [100 * en[en_map[c]]["prophage_cds_fraction_mean"] for c in cats]
bg = [100 * en[en_map[c]]["chromosome_cds_fraction_mean"] for c in cats]
ax.bar(x - 0.19, bg, width=0.36, color="#bdbdbd", label="Chromosome")
ax.bar(x + 0.19, pp, width=0.36, color="#c0392b", label="Prophage")
for i, c in enumerate(cats):
    q = en[en_map[c]]["q_value_bh"]
    star = "***" if q < 0.001 else ("**" if q < 0.01 else ("*" if q < 0.05 else "ns"))
    ax.text(i, max(pp[i], bg[i]) + 0.15, star, ha="center", fontsize=8)
ax.set_xticks(x)
ax.set_xticklabels(["Defence", "Virulence", "AMR", "MGE"])
ax.set_ylabel("% of CDS")
ax.set_ylim(0, max(max(pp), max(bg)) * 1.3)
ax.legend(frameon=False, loc="upper right", fontsize=8)
tag(ax, "b")

ax = fig.add_subplot(gs[1, :])
sys_path = os.path.join(RESULT, "defensefinder_systems.csv")
systems = pd.read_csv(sys_path)
systems["system_label"] = systems.apply(
    lambda r: "CRISPR-Cas" if r["model_family"] == "Cas" else r["system_name"], axis=1)
tab = systems.groupby("system_label").agg(complete=("complete", "sum"),
                                          partial=("partial", "sum"))
tab["total"] = tab["complete"] + tab["partial"]
# keep every type with at least one complete or partial system (13 complete types + 1 partial type)
tab = tab[tab["total"] > 0].sort_values(["complete", "total"], ascending=False)
ypos = np.arange(len(tab))[::-1]
ax.barh(ypos, tab["complete"].values, color="#c0392b", height=0.62, label="Complete")
ax.barh(ypos, tab["partial"].values, left=tab["complete"].values,
        color="#e8a49c", height=0.62, label="Partial")
ax.set_yticks(ypos)
ax.set_yticklabels(tab.index)
ax.set_xlabel(f"Defence systems ({int(tab['complete'].sum())} complete, {int(tab['partial'].sum())} partial)")
ax.legend(frameon=False, loc="lower right", fontsize=8)
ax.set_xlim(0, max(tab["total"].max() * 1.25, 5))
for y, (c, p) in zip(ypos, zip(tab["complete"].values, tab["partial"].values)):
    ax.text(c + p + 0.15, y, f"{int(c)}/{int(c + p)}", va="center", fontsize=8)
tag(ax, "c")
save(fig, "fig6_cargo_defence")

# ---------- Fig S1 ----------
fig, axes = plt.subplots(1, 2, figsize=(W, 2.75))
ax = axes[0]
for c in clades:
    sub = genomes[genomes["mlst_clade"] == c]
    ax.scatter(sub["genome_size"] / 1e6, sub["prophage_percent"], s=8, color=CLADE_COLORS[c],
               edgecolor="white", linewidth=0.2, label=f"Clade {c}", alpha=0.8)
rho, p = stats.spearmanr(genomes["genome_size"], genomes["prophage_percent"])
ax.text(0.03, 0.95, f"Spearman $\\rho$ = {rho:.2f}, {fmt_p(p)}", transform=ax.transAxes, fontsize=8, va="top")
ax.set_xlabel("Genome size (Mb)")
ax.set_ylabel("Prophage DNA (% of genome)")
ax.legend(frameon=False, loc="lower right", fontsize=8)
tag(ax, "a")

ax = axes[1]
for c in clades:
    sub = genomes[genomes["mlst_clade"] == c]
    ax.scatter(sub["genome_gc_percent"], sub["n_prophages"], s=8, color=CLADE_COLORS[c],
               edgecolor="white", linewidth=0.2, alpha=0.8)
rho2, p2 = stats.spearmanr(genomes["genome_gc_percent"], genomes["n_prophages"])
ax.text(0.97, 0.95, f"Spearman $\\rho$ = {rho2:.2f}, {fmt_p(p2)}", transform=ax.transAxes, fontsize=8,
        va="top", ha="right")
ax.set_xlabel("Genome GC (%)")
ax.set_ylabel("Prophage regions per genome")
tag(ax, "b")
fig.tight_layout()
save(fig, "figS1_genome_correlates")

# ---------- Fig S2 ----------
fig, axes = plt.subplots(1, 2, figsize=(W, 3.2))
ax = axes[0]
for conf in ["high", "medium", "low"]:
    sub = details[details["confidence"] == conf]
    ax.scatter(sub["best_qcov"] * 100, sub["best_pident"], s=7, color=CONF_COLORS[conf],
               edgecolor="white", linewidth=0.2, label=conf.capitalize(), alpha=0.75)
ax.axvline(50, color="#d62728", linestyle="--", linewidth=0.9)
ax.axhline(70, color="#d62728", linestyle="--", linewidth=0.9)
ax.text(52, 71.5, "typeable (≥70% id, ≥50% cov)", fontsize=8, color="#d62728", va="bottom",
        bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.5})
ax.set_xlabel("Best reference coverage (%)")
ax.set_ylabel("Best reference identity (%)")
ax.set_xlim(0, 105)
ax.set_ylim(60, 102)
ax.legend(frameon=False, markerscale=1.5, loc="lower right")
tag(ax, "a")

ax = axes[1]
typeable = details[details["phage_typeable"] == 1]["best_phage"].value_counts()
top = typeable.head(15)
ax.barh(range(len(top))[::-1], top.values, color="#2e75b6", height=0.65)
ax.set_yticks(range(len(top))[::-1])
ax.set_yticklabels(top.index)
ax.set_xlabel("Number of prophage regions")
tag(ax, "b")
fig.tight_layout()
save(fig, "figS2_phage_reference")

print("figures written")
