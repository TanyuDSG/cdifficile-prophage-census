#!/usr/bin/env python3
"""Publication figures for the all-241 prophage analysis."""
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
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from scipy import stats
from scipy.spatial.distance import pdist, squareform

BASE = "/home/bioinfo/cdprophage"
RESULT = os.path.join(BASE, "scaleup_241", "results")
OUT = os.path.join(BASE, "scaleup_241", "figures")
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.spines.top": False,
                     "axes.spines.right": False, "figure.dpi": 120})

CLADE_COLORS = {"1": "#1f77b4", "2": "#ff7f0e", "3": "#2ca02c", "4": "#d62728", "5": "#8c564b"}
CONF_COLORS = {"high": "#1f4e79", "medium": "#5b9bd5", "low": "#c9c9c9"}
PRESENT_COLOR = "#08306b"

genomes = pd.read_csv(os.path.join(RESULT, "genome_summary_analysis.csv"), dtype={"ST": str, "mlst_clade": str})
details = pd.read_csv(os.path.join(RESULT, "prophage_details.csv"), dtype={"ST": str, "mlst_clade": str})


def save(fig, name):
    fig.savefig(os.path.join(OUT, name + ".png"), dpi=300, bbox_inches="tight")
    fig.savefig(os.path.join(OUT, name + ".pdf"), bbox_inches="tight")
    plt.close(fig)


# ---------- Fig 1 ----------
fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.9), gridspec_kw={"width_ratios": [1.35, 1]})
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
    ax.scatter(x, vals, s=11, alpha=0.65, color=CLADE_COLORS[c], edgecolor="white", linewidth=0.2, zorder=3)
ax.set_xticks(range(len(clades)))
ax.set_xticklabels([f"Clade {c}\n(n={len(genomes[genomes['mlst_clade']==c])})" for c in clades], fontsize=7.5)
ax.set_ylabel("Prophage regions per genome")
ax.set_title("A  Prophage carriage by MLST clade", loc="left", fontweight="bold")

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
ax.set_yticklabels(labels, fontsize=8)
ax.set_xlim(0, 110)
ax.set_xlabel("%% of genomes (n=%d)" % len(genomes))
ax.set_title("B  Genome-level detection rates", loc="left", fontweight="bold")
fig.tight_layout()
save(fig, "fig1_carriage_by_clade")

# ---------- Fig 2 ----------
fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.7))
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
ax.set_xticklabels([f"{int(x/1000)}" for x in np.arange(0, bins[-1] + 1, 10000)], fontsize=7)
ax.legend(frameon=False, title="Confidence", fontsize=7, title_fontsize=7)
ax.set_title("A  Size distribution", loc="left", fontweight="bold")

ax = axes[1]
for conf in ["high", "medium", "low"]:
    sub = details[details["confidence"] == conf]
    ax.scatter(sub["length_bp"] / 1000, sub["gc_percent"], s=12, color=CONF_COLORS[conf],
               edgecolor="white", linewidth=0.2, label=conf.capitalize(), alpha=0.75)
mean_gc = genomes["genome_gc_percent"].mean()
ax.axhline(mean_gc, color="#d62728", linestyle="--", linewidth=1,
           label=f"Mean host GC ({mean_gc:.1f}%)")
ax.set_xlabel("Prophage region length (kb)")
ax.set_ylabel("Prophage GC (%)")
ax.legend(frameon=False, fontsize=6.2, markerscale=1.5, loc="upper right")
ax.set_title("B  Length vs GC", loc="left", fontweight="bold")
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

fig, axes = plt.subplots(1, 3, figsize=(11.0, 3.5), gridspec_kw={"width_ratios": [1.1, 0.9, 1.2]})
ax = axes[0]
ax.bar(range(1, len(counts) + 1), counts.values, color="#2e75b6", width=0.9)
ax.set_xlabel("Cluster rank")
ax.set_ylabel("Number of genomes carrying cluster")
ax.set_title("A  Cluster prevalence", loc="left", fontweight="bold")

ax = axes[1]
order = ["1", "2-3", "4-10", ">10"]
vals = [sharing.get(k, 0) for k in order]
ax.bar(range(len(order)), vals, color=["#c9c9c9", "#9dc3e6", "#2e75b6", "#1f4e79"], width=0.65)
for i, v in enumerate(vals):
    ax.text(i, v + max(vals) * 0.02, str(v), ha="center", fontsize=8)
ax.set_xticks(range(len(order)))
ax.set_xticklabels(["1", "2-3", "4-10", ">10"])
ax.set_ylabel("Number of clusters")
ax.set_title("B  Sharing breadth", loc="left", fontweight="bold")

ax = axes[2]
acc = pd.read_csv(os.path.join(RESULT, "cluster_accumulation.csv"))
ax.plot(acc["n_genomes"], acc["clusters_mean"], color="#1f4e79", linewidth=1.8)
ax.fill_between(acc["n_genomes"], acc["clusters_p05"], acc["clusters_p95"], color="#9dc3e6", alpha=0.35)
ax.set_xlabel("Number of genomes sampled")
ax.set_ylabel("Cumulative clusters")
ax.set_title("C  Cluster accumulation", loc="left", fontweight="bold")
fig.tight_layout()
save(fig, "fig3_cluster_sharing")

# ---------- Fig 4 ----------
tree = Phylo.read(os.path.join(RESULT, "mash_nj_tree.nwk"), "newick")
tree_order = [t.name for t in tree.get_terminals() if t.name in set(genomes["assembly_accession"])]
missing = [a for a in genomes["assembly_accession"] if a not in tree_order]
tree_order += missing
meta = genomes.set_index("assembly_accession").loc[tree_order]
top = list(counts.index[:60])
mat = presence.loc[tree_order, top].values
clade_codes = np.array([int(c) for c in meta["mlst_clade"]])

fig = plt.figure(figsize=(7.6, 7.8))
gs = fig.add_gridspec(2, 2, width_ratios=[0.045, 1], height_ratios=[0.10, 1],
                      hspace=0.02, wspace=0.01)
ax_top = fig.add_subplot(gs[0, 1])
ax_side = fig.add_subplot(gs[1, 0])
ax = fig.add_subplot(gs[1, 1])

ax_top.bar(range(len(top)), counts.loc[top].values, color="#4292c6", width=0.85)
ax_top.set_xlim(-0.5, len(top) - 0.5)
ax_top.set_ylabel("Genomes", fontsize=6)
ax_top.tick_params(axis="x", labelbottom=False, length=0)
ax_top.tick_params(axis="y", labelsize=5.5)
for s in ("top", "right", "bottom"):
    ax_top.spines[s].set_visible(False)
ax_top.set_title("Presence of the 60 most prevalent prophage clusters across 241 genomes",
                 loc="left", fontweight="bold", fontsize=9)

ax_side.imshow(clade_codes.reshape(-1, 1), aspect="auto", interpolation="nearest",
               cmap=ListedColormap([CLADE_COLORS[str(i)] for i in range(1, 6)]), vmin=1, vmax=5)
ax_side.set_xticks([]); ax_side.set_yticks([])
ax_side.set_ylabel("Genomes ordered by Mash tree", fontsize=7.5)
for spine in ax_side.spines.values():
    spine.set_visible(False)

ax.imshow(mat, aspect="auto", interpolation="nearest",
          cmap=ListedColormap(["#f7f7f7", PRESENT_COLOR]), vmin=0, vmax=1)
ax.set_xticks(range(len(top)))
ax.set_xticklabels(top, rotation=90, fontsize=3.4)
ax.set_yticks([])
ax.set_xlabel("Prophage cluster (ordered by prevalence)", fontsize=7.5)
ax.tick_params(axis="x", length=1.5, pad=1)
for spine in ax.spines.values():
    spine.set_linewidth(0.4)

legend_elements = [Patch(facecolor=CLADE_COLORS[str(i)], label=f"Clade {i}") for i in range(1, 6)]
ax.legend(handles=legend_elements, loc="upper left", bbox_to_anchor=(1.012, 1.0),
          frameon=False, fontsize=6.5, title="MLST clade", title_fontsize=6.5)
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
fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.7), gridspec_kw={"width_ratios": [1.25, 1, 0.95]})
ax = axes[0]
for c in clades:
    sel = meta["mlst_clade"].astype(str) == c
    ax.scatter(coords[sel, 0], coords[sel, 1], s=9, color=CLADE_COLORS[c], alpha=0.7,
               edgecolor="white", linewidth=0.15, label=f"Clade {c} (n={sel.sum()})")
ax.set_xlabel(f"PCoA1 ({explained[0]*100:.1f}% of variance)")
ax.set_ylabel(f"PCoA2 ({explained[1]*100:.1f}% of variance)")
ax.legend(frameon=False, fontsize=5.5, markerscale=1.5)
ax.set_title("A  Prophage gene-content space", loc="left", fontweight="bold")

ax = axes[1]
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
ax.text(0.03, 0.95, f"mean {np.mean(within):.3f} vs {np.mean(between):.3f}\nMann-Whitney P = {pw:.1e}",
        transform=ax.transAxes, fontsize=6.2, va="top")
json.dump({"within_mean": round(float(np.mean(within)), 4),
           "between_mean": round(float(np.mean(between)), 4),
           "n_pairs_within": int(len(within)), "n_pairs_between": int(len(between)),
           "mannwhitney_p": float(pw)},
          open(os.path.join(RESULT, "gene_sharing_stats.json"), "w"), indent=2)
ax.set_title("B  Gene sharing by clade", loc="left", fontweight="bold")

ax = axes[2]
pn = json.load(open(os.path.join(RESULT, "protein_network_summary.json")))
cats = ["phage_core", "defence", "MGE", "virulence", "AMR"]
labels = {"phage_core": "Phage core", "defence": "Defence", "MGE": "MGE",
          "virulence": "Virulence", "AMR": "AMR"}
tot = pn["gene_families_by_category"]
cross = pn["cross_clade_gene_families"]
vals2 = [100 * cross[c] / tot[c] if tot[c] else 0 for c in cats]
colors = ["#1f4e79", "#c0392b", "#7f8c8d", "#e67e22", "#2e86c1"]
ax.bar(range(len(cats)), vals2, color=colors, width=0.65)
for i, (v, c) in enumerate(zip(vals2, cats)):
    ax.text(i, v + 1.5, f"{v:.0f}%\n({cross[c]}/{tot[c]})", ha="center", fontsize=5.8)
ax.set_xticks(range(len(cats)))
ax.set_xticklabels([labels[c] for c in cats], rotation=30, ha="right", fontsize=7)
ax.set_ylabel("% of gene families in >=2 clades")
ax.set_ylim(0, max(60, max(vals2) * 1.25))
ax.set_title("C  Cross-clade sharing", loc="left", fontweight="bold")
fig.tight_layout()
save(fig, "fig5_protein_network")

# ---------- Fig 6 ----------
cargo = pd.read_csv(os.path.join(RESULT, "prophage_cargo_summary.csv"), dtype={"ST": str, "mlst_clade": str})
cs = json.load(open(os.path.join(RESULT, "cargo_summary.json")))
fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.8), gridspec_kw={"width_ratios": [1, 1, 1.15]})
ax = axes[0]
cats = ["defence", "virulence", "AMR", "MGE"]
labels = ["Defence", "Virulence", "AMR", "MGE"]
colmap = {"defence": "has_defence", "virulence": "has_virulence", "AMR": "has_amr", "MGE": "has_mge"}
high_vals = [100 * cargo[cargo["confidence"] == "high"][colmap[c]].mean() for c in cats]
all_vals = [100 * cargo[colmap[c]].mean() for c in cats]
x = np.arange(len(cats))
ax.bar(x - 0.19, all_vals, width=0.36, color="#9dc3e6", label=f"All regions (n={len(cargo)})")
ax.bar(x + 0.19, high_vals, width=0.36, color="#1f4e79", label=f"High-confidence (n={(cargo['confidence']=='high').sum()})")
for i, (a, h) in enumerate(zip(all_vals, high_vals)):
    ax.text(i - 0.19, a + 1.2, f"{a:.0f}%", ha="center", fontsize=6.2)
    ax.text(i + 0.19, h + 1.2, f"{h:.0f}%", ha="center", fontsize=6.2)
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.set_ylabel("% of prophage regions")
ax.set_ylim(0, 100)
ax.legend(frameon=False, fontsize=6.2, loc="upper right")
ax.set_title("A  Cargo prevalence", loc="left", fontweight="bold")

ax = axes[1]
en = cs["enrichment_vs_chromosome"]
en_map = {"defence": "defence", "virulence": "virulence", "AMR": "amr", "MGE": "mge"}
x = np.arange(len(cats))
pp = [100 * en[en_map[c]]["prophage_cds_fraction_mean"] for c in cats]
bg = [100 * en[en_map[c]]["chromosome_cds_fraction_mean"] for c in cats]
ax.bar(x - 0.19, bg, width=0.36, color="#bdbdbd", label="Chromosome (prophage-excluded)")
ax.bar(x + 0.19, pp, width=0.36, color="#c0392b", label="Prophage regions")
for i, c in enumerate(cats):
    q = en[en_map[c]]["q_value_bh"]
    star = "***" if q < 0.001 else ("**" if q < 0.01 else ("*" if q < 0.05 else "ns"))
    ax.text(i, max(pp[i], bg[i]) + 0.08, star, ha="center", fontsize=8)
ax.set_xticks(x)
ax.set_xticklabels(["Defence", "Virulence", "AMR", "MGE"])
ax.set_ylabel("% of CDS")
ax.set_ylim(0, max(max(pp), max(bg)) * 1.25)
ax.legend(frameon=False, fontsize=5.8, loc="upper right")
ax.set_title("B  Enrichment vs chromosome", loc="left", fontweight="bold")

ax = axes[2]
sys_path = os.path.join(RESULT, "defensefinder_systems.csv")
if os.path.exists(sys_path):
    systems = pd.read_csv(sys_path)
    systems["system_label"] = systems.apply(
        lambda r: "CRISPR-Cas" if r["model_family"] == "Cas" else r["system_name"], axis=1)
    tab = systems.groupby("system_label").agg(complete=("complete", "sum"),
                                              partial=("partial", "sum"))
    tab["total"] = tab["complete"] + tab["partial"]
    tab = tab.sort_values("total", ascending=False).head(12)
    ypos = np.arange(len(tab))[::-1]
    ax.barh(ypos, tab["complete"].values, color="#c0392b", height=0.62, label="Complete")
    ax.barh(ypos, tab["partial"].values, left=tab["complete"].values,
            color="#e8a49c", height=0.62, label="Partial")
    ax.set_yticks(ypos)
    ax.set_yticklabels(tab.index, fontsize=6.2)
    ax.set_xlabel("Defence systems (regions)")
    ax.legend(frameon=False, fontsize=6, loc="lower right")
    ax.set_xlim(0, max(tab["total"].max() * 1.15, 5))
    for y, (c, p) in zip(ypos, zip(tab["complete"].values, tab["partial"].values)):
        if c + p > 0:
            ax.text(c + p + 0.5, y, f"{int(c)}/{int(c + p)}", va="center", fontsize=5.6)
else:
    dc = cs["defence_category_counts"]
    dc = dict(sorted(dc.items(), key=lambda kv: kv[1], reverse=True))
    ax.barh(range(len(dc))[::-1], list(dc.values()), color="#c0392b", height=0.62)
    ax.set_yticks(range(len(dc))[::-1])
    ax.set_yticklabels(list(dc.keys()), fontsize=7)
ax.set_title("C  Defence systems reconstructed", loc="left", fontweight="bold")
fig.tight_layout()
save(fig, "fig6_cargo_defence")

# ---------- Fig S1 ----------
fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.6))
ax = axes[0]
for c in clades:
    sub = genomes[genomes["mlst_clade"] == c]
    ax.scatter(sub["genome_size"] / 1e6, sub["prophage_percent"], s=12, color=CLADE_COLORS[c],
               edgecolor="white", linewidth=0.2, label=f"Clade {c}", alpha=0.8)
rho, p = stats.spearmanr(genomes["genome_size"], genomes["prophage_percent"])
ax.text(0.03, 0.95, f"Spearman rho = {rho:.2f}, P = {p:.2g}", transform=ax.transAxes, fontsize=7, va="top")
ax.set_xlabel("Genome size (Mb)")
ax.set_ylabel("Prophage DNA (% of genome)")
ax.legend(frameon=False, fontsize=6, loc="lower right")
ax.set_title("A  Genome size vs prophage burden", loc="left", fontweight="bold")

ax = axes[1]
for c in clades:
    sub = genomes[genomes["mlst_clade"] == c]
    ax.scatter(sub["genome_gc_percent"], sub["n_prophages"], s=12, color=CLADE_COLORS[c],
               edgecolor="white", linewidth=0.2, alpha=0.8)
rho2, p2 = stats.spearmanr(genomes["genome_gc_percent"], genomes["n_prophages"])
ax.text(0.03, 0.95, f"Spearman rho = {rho2:.2f}, P = {p2:.2g}", transform=ax.transAxes, fontsize=7, va="top")
ax.set_xlabel("Genome GC (%)")
ax.set_ylabel("Prophage regions per genome")
ax.set_title("B  GC vs prophage count", loc="left", fontweight="bold")
fig.tight_layout()
save(fig, "figS1_genome_correlates")

# ---------- Fig S2 ----------
fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.6))
ax = axes[0]
for conf in ["high", "medium", "low"]:
    sub = details[details["confidence"] == conf]
    ax.scatter(sub["best_qcov"] * 100, sub["best_pident"], s=10, color=CONF_COLORS[conf],
               edgecolor="white", linewidth=0.2, label=conf.capitalize(), alpha=0.75)
ax.axvline(50, color="#d62728", linestyle="--", linewidth=0.9)
ax.axhline(70, color="#d62728", linestyle="--", linewidth=0.9)
ax.text(51, 71, "typeable\n(>=70% id, >=50% cov)", fontsize=6.2, color="#d62728", va="bottom")
ax.set_xlabel("Best reference coverage (%)")
ax.set_ylabel("Best reference identity (%)")
ax.set_xlim(0, 105)
ax.set_ylim(60, 102)
ax.legend(frameon=False, fontsize=6.5, markerscale=1.5)
ax.set_title("A  Similarity to complete phages", loc="left", fontweight="bold")

ax = axes[1]
typeable = details[details["phage_typeable"] == 1]["best_phage"].value_counts()
top = typeable.head(15)
ax.barh(range(len(top))[::-1], top.values, color="#2e75b6", height=0.65)
ax.set_yticks(range(len(top))[::-1])
ax.set_yticklabels(top.index, fontsize=6)
ax.set_xlabel("Number of prophage regions")
ax.set_title("B  Closest reference phage", loc="left", fontweight="bold")
fig.tight_layout()
save(fig, "figS2_phage_reference")

print("figures written")
