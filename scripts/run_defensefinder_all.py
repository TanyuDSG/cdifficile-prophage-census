#!/usr/bin/env python3
"""Run DefenseFinder/MacSyFinder on all prophage proteins and merge results.

Three model sets are run independently in unordered mode:
  defense-finder-models/DefenseFinder
  defense-finder-models/RM
  CasFinder

MacSyFinder 2.x writes gene-level candidate systems to all_systems.tsv in
unordered mode. Systems with sys_wholeness == 1.0 are treated as complete;
0.5-0.99 are reported as partial.
"""
import json
import os
import subprocess
from collections import Counter, defaultdict

import pandas as pd

BASE = "/home/bioinfo/cdprophage"
RESULT = os.path.join(BASE, "scaleup_241", "results")
OUT = os.path.join(RESULT, "defensefinder")
PROTEINS = os.path.join(RESULT, "prophage_proteins.fasta")
LOG = os.path.join(BASE, "scaleup_241", "logs")
os.makedirs(OUT, exist_ok=True)
os.makedirs(LOG, exist_ok=True)

ENV = os.environ.copy()
ENV["PATH"] = os.pathsep.join([
    os.path.join(BASE, "tools", "hmmer", "usr", "bin"),
    os.path.expanduser("~/.local/bin"),
    "/usr/local/sbin", "/usr/local/bin", "/usr/sbin", "/usr/bin", "/sbin", "/bin",
])
MACSY = os.path.expanduser("~/.local/bin/macsyfinder")

MODELS = [
    ("DefenseFinder", ["defense-finder-models/DefenseFinder", "all"],
     ["--coverage-profile", "0.4", "--exchangeable-weight", "1"]),
    ("RM", ["defense-finder-models/RM", "all"],
     ["--coverage-profile", "0.4", "--exchangeable-weight", "1"]),
    ("Cas", ["CasFinder", "all"], []),
]


def run_models():
    if not os.path.exists(PROTEINS):
        raise SystemExit("protein fasta not found: " + PROTEINS)
    for name, models, extra in MODELS:
        outdir = os.path.join(OUT, name)
        tsv = os.path.join(outdir, "all_systems.tsv")
        if os.path.exists(tsv):
            continue
        cmd = [MACSY, "--db-type", "unordered", "--sequence-db", PROTEINS,
               "--models"] + models + ["--out-dir", outdir, "--w", "8"] + extra
        with open(os.path.join(LOG, f"macsyfinder_{name}.log"), "w") as fh:
            subprocess.run(cmd, check=True, stdout=fh, stderr=subprocess.STDOUT, env=ENV)


def parse():
    frames = []
    for name, _, _ in MODELS:
        path = os.path.join(OUT, name, "all_systems.tsv")
        if not os.path.exists(path):
            continue
        try:
            df = pd.read_csv(path, sep="\t", comment="#")
        except pd.errors.EmptyDataError:
            continue
        if df.empty:
            continue
        df["model_family"] = name
        frames.append(df)
    if not frames:
        return pd.DataFrame(), pd.DataFrame()
    genes = pd.concat(frames, ignore_index=True)
    genes["prophage_id"] = genes["hit_id"].astype(str).str.split("|").str[0]
    genes["system_name"] = genes["model_fqn"].astype(str).str.split("/").str[-2]
    genes["subtype"] = genes["model_fqn"].astype(str).str.split("/").str[-1]
    systems = genes.groupby("sys_id", as_index=False).agg(
        model_family=("model_family", "first"),
        system_name=("system_name", "first"),
        subtype=("subtype", "first"),
        wholeness=("sys_wholeness", "first"),
        n_genes=("hit_id", "nunique"),
        prophage_id=("prophage_id", "first"),
        genes=("hit_id", lambda s: ";".join(sorted(set(s)))),
    )
    systems["complete"] = (systems["wholeness"] >= 0.999).astype(int)
    systems["partial"] = ((systems["wholeness"] >= 0.5) & (systems["wholeness"] < 0.999)).astype(int)
    return genes, systems


def main():
    run_models()
    genes, systems = parse()
    if genes.empty:
        raise SystemExit("no DefenseFinder hits parsed")
    genes.to_csv(os.path.join(RESULT, "defensefinder_genes.csv"), index=False)
    systems.to_csv(os.path.join(RESULT, "defensefinder_systems.csv"), index=False)

    cargo = pd.read_csv(os.path.join(RESULT, "prophage_cargo_summary.csv"),
                        dtype={"ST": str, "mlst_clade": str})
    per = defaultdict(lambda: {"complete": 0, "partial": 0, "types": set()})
    for _, r in systems.iterrows():
        d = per[r["prophage_id"]]
        if r["complete"] and r["model_family"] != "Cas":
            d["complete"] += 1
            d["types"].add(r["system_name"])
        elif r["complete"] and r["model_family"] == "Cas":
            d["complete"] += 1
            d["types"].add("CRISPR-Cas")
        elif r["partial"]:
            d["partial"] += 1
    cargo["n_df_systems"] = cargo["prophage_id"].map(lambda p: per[p]["complete"]).fillna(0).astype(int)
    cargo["n_df_partial"] = cargo["prophage_id"].map(lambda p: per[p]["partial"]).fillna(0).astype(int)
    cargo["df_system_types"] = cargo["prophage_id"].map(
        lambda p: ";".join(sorted(per[p]["types"]))).fillna("")
    cargo["has_df_system"] = (cargo["n_df_systems"] > 0).astype(int)
    cargo["has_any_df_hit"] = ((cargo["n_df_systems"] + cargo["n_df_partial"]) > 0).astype(int)
    cargo.to_csv(os.path.join(RESULT, "prophage_cargo_summary.csv"), index=False)

    complete = systems[systems["complete"] == 1]
    summary = {
        "proteins_searched": int(genes["hit_id"].nunique()),
        "systems_total": int(len(systems)),
        "systems_complete": int(len(complete)),
        "systems_partial": int(systems["partial"].sum()),
        "prophages_with_complete_system": int((cargo["n_df_systems"] > 0).sum()),
        "prophages_with_any_hit": int((cargo["has_any_df_hit"] > 0).sum()),
        "prophage_complete_system_rate": round(float((cargo["n_df_systems"] > 0).mean()), 4),
        "family_counts_complete": dict(Counter(complete["model_family"])),
        "system_type_counts_complete": dict(Counter(complete["system_name"])),
        "system_type_counts_all": dict(Counter(systems["system_name"])),
    }
    with open(os.path.join(RESULT, "defensefinder_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)

    cpath = os.path.join(RESULT, "cargo_summary.json")
    if os.path.exists(cpath):
        with open(cpath) as fh:
            cs = json.load(fh)
        cs["defensefinder"] = summary
        cs["defence_any_source_rate_all"] = round(float(
            ((cargo["has_defence"] == 1) | (cargo["has_any_df_hit"] == 1)).mean()), 4)
        cs["defence_any_source_rate_high"] = round(float(
            ((cargo[cargo["confidence"] == "high"]["has_defence"] == 1) |
             (cargo[cargo["confidence"] == "high"]["has_any_df_hit"] == 1)).mean()), 4)
        with open(cpath, "w") as fh:
            json.dump(cs, fh, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
