#!/usr/bin/env python3
"""
argo_amr_summary.py

Comprehensive summary of ARGO AMR results across all samples from:
  ./Argo_results

Output directory: ./argo_amr_summary_output/

Outputs:
  per_sample_stats.csv       – richness, diversity, plasmid fraction per sample
  gene_frequency.csv         – prevalence + classification of each gene
  drug_class_profile.csv     – sample × drug_class copy-number matrix
  presence_absence.csv       – binary sample × gene matrix
  jaccard_similarity.csv     – pairwise Jaccard similarity between samples
  top_lineages.csv           – top ARG-carrying taxa per drug class
  01_gene_prevalence.png     – top genes ranked by % of samples containing them
  02_drug_class_heatmap.png  – drug class × sample log-abundance heatmap (clustered)
  03_jaccard_heatmap.png     – sample × sample Jaccard similarity heatmap (clustered)
  04_plasmid_burden.png      – plasmid vs chromosomal ARG fraction per sample
  05_sample_richness.png     – ARG gene richness and drug class count per sample
"""

import os, sys, json, csv, re
from collections import defaultdict, Counter
import numpy as np

try:
    import pandas as pd
    HAS_PD = True
except ImportError:
    HAS_PD = False
    print("[WARN] pandas not found; falling back to csv writer")

try:
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    HAS_MPL = True
except ImportError:
    HAS_MPL = False
    print("[WARN] matplotlib not found; plots will be skipped")

try:
    from scipy.cluster.hierarchy import linkage, leaves_list
    from scipy.spatial.distance import squareform
    HAS_SCI = True
except ImportError:
    HAS_SCI = False
    print("[WARN] scipy not found; plots will not be clustered")

# ── Directories & thresholds ──────────────────────────────────────────────────
ARGO_DIRS = [
    "./Argo_results"
]
OUT_DIR = "argo_amr_summary_output"

CORE_THR   = 0.80   # gene in >= 80% of samples
COMMON_THR = 0.50
RARE_THR   = 0.20   # gene in < 20% = rare
# ─────────────────────────────────────────────────────────────────────────────


# ── Sample metadata ────────────────────────────────────────────────────────────
def parse_sample_name(name):
    """Extract site and timepoint from names like Gentry_26_MR1015E."""
    m = re.search(r'_([A-Z]{2,}[\d]+)([A-Z])$', name)
    if m:
        return m.group(1), m.group(2)          # site, timepoint
    m = re.search(r'barcode(\d+)', name)
    if m:
        return f"barcode{m.group(1)}", None
    return name, None


# ── File discovery ─────────────────────────────────────────────────────────────
def find_samples(argo_dirs):
    samples = []
    for base in argo_dirs:
        if not os.path.isdir(base):
            print(f"[WARN] Not found: {base}")
            continue
        for sub in sorted(os.listdir(base)):
            dp = os.path.join(base, sub)
            if not os.path.isdir(dp):
                continue
            files = os.listdir(dp)
            tsvs  = [f for f in files if f.endswith(".sarg.tsv")]
            jsons = [f for f in files if f.endswith(".sarg.json")]
            if not tsvs or not jsons:
                continue
            name = tsvs[0].replace(".sarg.tsv", "").replace("chopped_", "")
            samples.append({
                "name":      name,
                "tsv_path":  os.path.join(dp, tsvs[0]),
                "json_path": os.path.join(dp, jsons[0]),
            })
    return samples


# ── Data loading ───────────────────────────────────────────────────────────────
def load_tsv(path):
    """
    Load .sarg.tsv → list of dicts.
    Columns: lineage, type, subtype, carrier, copy, genome, abundance
    """
    rows = []
    try:
        with open(path) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                try:
                    row["copy"]      = float(row.get("copy",      0) or 0)
                    row["abundance"] = float(row.get("abundance", 0) or 0)
                except (ValueError, TypeError):
                    row["copy"] = row["abundance"] = 0.0
                rows.append(row)
    except Exception as e:
        print(f"    [WARN] Could not read {path}: {e}")
    return rows


def load_json_read_counts(path):
    """
    Load .sarg.json → per-gene read counts + plasmid read counts.
    Returns gene_reads: Counter, plasmid_reads: Counter
    """
    gene_reads    = Counter()
    plasmid_reads = Counter()
    try:
        with open(path) as fh:
            data = json.load(fh)
        if not isinstance(data, dict):
            return gene_reads, plasmid_reads
        for info in data.values():
            if info.get("remark") != "ARG-containing":
                continue
            is_plasmid = bool(info.get("plasmid", False))
            for hit_str in info.get("hit", []):
                parts = hit_str.split("|")
                gene  = parts[1].strip() if len(parts) > 1 else hit_str
                gene_reads[gene]    += 1
                if is_plasmid:
                    plasmid_reads[gene] += 1
    except Exception as e:
        print(f"    [WARN] Could not read {path}: {e}")
    return gene_reads, plasmid_reads


# ── Per-sample statistics ──────────────────────────────────────────────────────
def per_sample_stats(samples_data):
    rows = []
    for s in samples_data:
        tsv  = s["tsv"]
        site, tp = parse_sample_name(s["name"])
        if not tsv:
            continue

        genes        = {r["subtype"] for r in tsv}
        drug_classes = {r["type"]    for r in tsv}
        total_copy   = sum(r["copy"] for r in tsv)
        plas_copy    = sum(r["copy"] for r in tsv
                          if r.get("carrier", "").lower() == "plasmid")

        # Shannon diversity of drug classes (by copy number)
        class_copies = defaultdict(float)
        for r in tsv:
            class_copies[r["type"]] += r["copy"]
        total = sum(class_copies.values())
        shannon = 0.0
        for c in class_copies.values():
            p = c / total if total > 0 else 0
            if p > 0:
                shannon -= p * np.log2(p)

        gr  = s["gene_reads"]
        pr  = s["plasmid_reads"]
        tr  = sum(gr.values())
        plr = sum(pr.values())

        rows.append({
            "sample":                  s["name"],
            "site":                    site,
            "timepoint":               tp or "",
            "n_unique_genes":          len(genes),
            "n_drug_classes":          len(drug_classes),
            "total_arg_copies":        round(total_copy, 2),
            "plasmid_arg_copies":      round(plas_copy,  2),
            "plasmid_frac_copies":     round(plas_copy / total_copy, 4) if total_copy > 0 else 0,
            "total_arg_reads":         tr,
            "plasmid_arg_reads":       plr,
            "plasmid_frac_reads":      round(plr / tr, 4) if tr > 0 else 0,
            "shannon_diversity_bits":  round(shannon, 4),
        })
    return rows


# ── Gene frequency ─────────────────────────────────────────────────────────────
def gene_frequency_table(samples_data, n_samples):
    gene_info = defaultdict(lambda: {
        "drug_class": "", "n_samples": 0,
        "total_copies": 0.0, "total_reads": 0,
        "plasmid_reads": 0, "n_plasmid_samples": 0
    })

    for s in samples_data:
        seen = set()
        for r in s["tsv"]:
            g = r["subtype"]
            gene_info[g]["drug_class"]   = r["type"]
            gene_info[g]["total_copies"] += r["copy"]
            if g not in seen:
                gene_info[g]["n_samples"] += 1
                seen.add(g)
        for gene, cnt in s["gene_reads"].items():
            gene_info[gene]["total_reads"]   += cnt
        for gene, cnt in s["plasmid_reads"].items():
            gene_info[gene]["plasmid_reads"] += cnt

    # Plasmid sample count: samples where this gene's plasmid reads > 0
    for s in samples_data:
        for gene in s["plasmid_reads"]:
            gene_info[gene]["n_plasmid_samples"] += 1

    rows = []
    for gene, info in sorted(gene_info.items(),
                              key=lambda x: x[1]["n_samples"], reverse=True):
        frac = info["n_samples"] / n_samples if n_samples > 0 else 0
        if frac >= CORE_THR:
            category = "core"
        elif frac >= COMMON_THR:
            category = "common"
        elif frac >= RARE_THR:
            category = "intermediate"
        else:
            category = "rare"

        plasmid_frac = (info["plasmid_reads"] / info["total_reads"]
                        if info["total_reads"] > 0 else 0)
        rows.append({
            "gene":              gene,
            "drug_class":        info["drug_class"],
            "n_samples":         info["n_samples"],
            "pct_samples":       round(100 * frac, 1),
            "prevalence_class":  category,
            "total_copies":      round(info["total_copies"], 2),
            "total_reads":       info["total_reads"],
            "plasmid_frac_reads":round(plasmid_frac, 4),
            "n_plasmid_samples": info["n_plasmid_samples"],
        })
    return rows


# ── Presence/absence matrix ────────────────────────────────────────────────────
def build_pa_matrix(samples_data):
    all_genes = sorted({r["subtype"]
                        for s in samples_data for r in s["tsv"]})
    sample_names = [s["name"] for s in samples_data]
    gene_index   = {g: i for i, g in enumerate(all_genes)}

    mat = np.zeros((len(sample_names), len(all_genes)), dtype=np.int8)
    for si, s in enumerate(samples_data):
        for r in s["tsv"]:
            mat[si, gene_index[r["subtype"]]] = 1
    return mat, sample_names, all_genes


# ── Jaccard matrix ─────────────────────────────────────────────────────────────
def jaccard_matrix(mat):
    n = mat.shape[0]
    J = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            inter = int(np.sum(mat[i] & mat[j]))
            union = int(np.sum(mat[i] | mat[j]))
            val   = inter / union if union > 0 else 0.0
            J[i, j] = J[j, i] = val
    return J


# ── Drug class profile matrix ──────────────────────────────────────────────────
def drug_class_matrix(samples_data):
    all_classes  = sorted({r["type"] for s in samples_data for r in s["tsv"]})
    sample_names = [s["name"] for s in samples_data]
    ci = {c: i for i, c in enumerate(all_classes)}

    mat = np.zeros((len(sample_names), len(all_classes)))
    for si, s in enumerate(samples_data):
        for r in s["tsv"]:
            mat[si, ci[r["type"]]] += r["copy"]
    return mat, sample_names, all_classes


# ── Top lineages ───────────────────────────────────────────────────────────────
def top_lineages(samples_data, top_n=10):
    lineage_class = defaultdict(lambda: defaultdict(float))
    for s in samples_data:
        for r in s["tsv"]:
            lineage_class[r["type"]][r["lineage"]] += r["copy"]

    rows = []
    for dc, lin_dict in sorted(lineage_class.items()):
        for lin, copies in sorted(lin_dict.items(), key=lambda x: x[1], reverse=True)[:top_n]:
            rows.append({"drug_class": dc, "lineage": lin, "total_copies": round(copies, 2)})
    return rows


# ── Hierarchical clustering order ─────────────────────────────────────────────
def cluster_order(mat):
    """Return sample order after average-linkage clustering of PA matrix."""
    if not HAS_SCI or mat.shape[0] < 3:
        return list(range(mat.shape[0]))
    dist = np.zeros((mat.shape[0], mat.shape[0]))
    for i in range(mat.shape[0]):
        for j in range(i, mat.shape[0]):
            union = int(np.sum(mat[i] | mat[j]))
            inter = int(np.sum(mat[i] & mat[j]))
            d = 1 - (inter / union if union > 0 else 0)
            dist[i, j] = dist[j, i] = d
    Z = linkage(squareform(dist), method="average")
    return list(leaves_list(Z))


# ── Plots ──────────────────────────────────────────────────────────────────────
DRUG_CLASS_COLORS = {
    "multidrug":    "#4e79a7",
    "tetracycline": "#f28e2b",
    "bacitracin":   "#76b7b2",
    "macrolide-lincosamide-streptogramin": "#59a14f",
    "beta-lactam":  "#e15759",
    "aminoglycoside":"#af7aa1",
    "sulfonamide":  "#ff9da7",
    "biocide":      "#9c755f",
    "colistin":     "#bab0ac",
    "defensin":     "#edc948",
    "rifamycin":    "#b07aa1",
    "phenicol":     "#d4a6c8",
    "glycopeptide": "#86bcb6",
    "quinolone":    "#f1ce63",
    "trimethoprim": "#499894",
}

def dc_color(dc):
    return DRUG_CLASS_COLORS.get(dc.lower(),
           DRUG_CLASS_COLORS.get(dc.lower().split()[0], "#aaaaaa"))


def plot_gene_prevalence(gene_freq_rows, n_samples, path, top_n=40):
    if not HAS_MPL: return
    rows = gene_freq_rows[:top_n]
    genes  = [r["gene"]        for r in rows]
    pcts   = [r["pct_samples"] for r in rows]
    dcs    = [r["drug_class"]  for r in rows]
    colors = [dc_color(d) for d in dcs]

    fig, ax = plt.subplots(figsize=(14, max(6, len(genes)*0.28)))
    ax.barh(range(len(genes)), pcts, color=colors, edgecolor="white", linewidth=0.4)
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontsize=8)
    ax.invert_yaxis()
    ax.axvline(80, color="red",    ls="--", lw=0.8, alpha=0.6, label="Core 80%")
    ax.axvline(50, color="orange", ls="--", lw=0.8, alpha=0.6, label="Common 50%")
    ax.axvline(20, color="gray",   ls="--", lw=0.8, alpha=0.4, label="Rare 20%")
    ax.set_xlabel("% of samples")
    ax.set_title(f"ARG Gene Prevalence — top {top_n} of {len(gene_freq_rows)} genes "
                 f"(n={n_samples} samples)", fontsize=11)
    ax.legend(fontsize=8, loc="lower right")

    # Drug class legend patches
    from matplotlib.patches import Patch
    seen_dc = {}
    for dc, col in zip(dcs, colors):
        seen_dc[dc] = col
    legend_els = [Patch(facecolor=v, label=k) for k, v in sorted(seen_dc.items())]
    ax.legend(handles=legend_els, fontsize=7, bbox_to_anchor=(1.01, 1),
              loc="upper left", title="Drug class")
    plt.tight_layout()
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {path}")


def plot_drug_class_heatmap(dc_mat, samples, drug_classes, order, path):
    if not HAS_MPL: return
    mat    = np.log1p(dc_mat[order])
    labels = [samples[i] for i in order]

    fig, ax = plt.subplots(figsize=(max(10, len(drug_classes)*0.85),
                                    max(8,  len(samples)*0.32)))
    im = ax.imshow(mat, aspect="auto", cmap="YlOrRd", interpolation="nearest")
    ax.set_xticks(range(len(drug_classes)))
    ax.set_xticklabels(drug_classes, rotation=45, ha="right", fontsize=9)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_title("ARG Drug-Class Abundance per Sample  (log₁₊ copy number, "
                 "hierarchically clustered)", fontsize=11)
    plt.colorbar(im, ax=ax, shrink=0.6, label="log₁₊(copies)")
    plt.tight_layout()
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {path}")


def plot_jaccard_heatmap(J, samples, order, path):
    if not HAS_MPL: return
    n   = len(samples)
    mat = J[np.ix_(order, order)]
    lbl = [samples[i] for i in order]

    fig, ax = plt.subplots(figsize=(max(10, n*0.38), max(8, n*0.38)))
    im = ax.imshow(mat, cmap="Blues", vmin=0, vmax=1,
                   aspect="auto", interpolation="nearest")
    ax.set_xticks(range(n)); ax.set_yticks(range(n))
    ax.set_xticklabels(lbl, rotation=90, fontsize=7)
    ax.set_yticklabels(lbl, fontsize=7)
    ax.set_title("Pairwise Jaccard Similarity of ARG Gene Sets  "
                 "(average-linkage cluster order)", fontsize=11)
    plt.colorbar(im, ax=ax, shrink=0.6, label="Jaccard similarity")
    plt.tight_layout()
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {path}")


def plot_plasmid_burden(stats_rows, path):
    if not HAS_MPL: return
    rows = sorted(stats_rows, key=lambda x: x["plasmid_frac_copies"], reverse=True)
    labels = [r["sample"] for r in rows]
    pf     = [r["plasmid_frac_copies"]         for r in rows]
    cf     = [1 - r["plasmid_frac_copies"]     for r in rows]
    x = np.arange(len(labels))

    fig, ax = plt.subplots(figsize=(max(14, len(labels)*0.45), 5))
    ax.bar(x, cf, label="Chromosomal", color="steelblue", alpha=0.85)
    ax.bar(x, pf, bottom=cf, label="Plasmid", color="tomato", alpha=0.85)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=90, fontsize=7)
    ax.set_ylabel("Fraction of ARG copies")
    ax.set_title("Plasmid vs Chromosomal ARG Fraction per Sample", fontsize=11)
    ax.legend()
    plt.tight_layout()
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {path}")


def plot_richness(stats_rows, order_names, path):
    if not HAS_MPL: return
    lookup = {r["sample"]: r for r in stats_rows}
    rows   = [lookup[n] for n in order_names if n in lookup]
    labels = [r["sample"]        for r in rows]
    genes  = [r["n_unique_genes"]  for r in rows]
    dcs    = [r["n_drug_classes"]  for r in rows]
    x = np.arange(len(labels))

    fig, ax1 = plt.subplots(figsize=(max(14, len(labels)*0.45), 5))
    ax2 = ax1.twinx()
    ax1.bar(x - 0.2, genes, 0.4, label="Unique genes",   color="steelblue", alpha=0.85)
    ax2.bar(x + 0.2, dcs,   0.4, label="Drug classes",   color="darkorange",alpha=0.85)
    ax1.set_xticks(x)
    ax1.set_xticklabels(labels, rotation=90, fontsize=7)
    ax1.set_ylabel("Unique ARG genes",  color="steelblue")
    ax2.set_ylabel("Drug classes",      color="darkorange")
    ax1.set_title("ARG Richness per Sample  (clustered order)", fontsize=11)
    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1+h2, l1+l2, fontsize=9)
    plt.tight_layout()
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {path}")


# ── CSV helpers ────────────────────────────────────────────────────────────────
def write_csv(rows, path, fieldnames=None):
    if not rows:
        return
    if fieldnames is None:
        fieldnames = list(rows[0].keys())
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)
    print(f"  Saved: {path}")


def write_matrix_csv(mat, row_labels, col_labels, path):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["sample"] + col_labels)
        for i, lbl in enumerate(row_labels):
            w.writerow([lbl] + list(mat[i]))
    print(f"  Saved: {path}")


# ── Main ───────────────────────────────────────────────────────────────────────
def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    # ── 1. Discover samples ────────────────────────────────────────────────────
    print("\n=== Discovering samples ===")
    sample_files = find_samples(ARGO_DIRS)
    if not sample_files:
        sys.exit("No samples found. Check ARGO_DIRS paths.")
    print(f"Found {len(sample_files)} samples\n")

    # ── 2. Load data ───────────────────────────────────────────────────────────
    print("=== Loading data ===")
    samples_data = []
    for sf in sample_files:
        tsv  = load_tsv(sf["tsv_path"])
        gr, pr = load_json_read_counts(sf["json_path"])
        samples_data.append({
            "name":          sf["name"],
            "tsv":           tsv,
            "gene_reads":    gr,
            "plasmid_reads": pr,
        })
        print(f"  {sf['name']:50s}  {len(tsv):4d} TSV rows  "
              f"{sum(gr.values()):5d} ARG reads")

    n_samples = len(samples_data)

    # ── 3. Per-sample stats ────────────────────────────────────────────────────
    print("\n=== Computing per-sample statistics ===")
    stats_rows = per_sample_stats(samples_data)
    write_csv(stats_rows, os.path.join(OUT_DIR, "per_sample_stats.csv"))

    # ── 4. Gene frequency ──────────────────────────────────────────────────────
    print("\n=== Computing gene frequency table ===")
    gf_rows = gene_frequency_table(samples_data, n_samples)
    write_csv(gf_rows, os.path.join(OUT_DIR, "gene_frequency.csv"))

    n_core   = sum(1 for r in gf_rows if r["prevalence_class"] == "core")
    n_common = sum(1 for r in gf_rows if r["prevalence_class"] == "common")
    n_inter  = sum(1 for r in gf_rows if r["prevalence_class"] == "intermediate")
    n_rare   = sum(1 for r in gf_rows if r["prevalence_class"] == "rare")
    print(f"  Total unique genes : {len(gf_rows)}")
    print(f"  Core   (>= {int(CORE_THR*100)}%): {n_core}")
    print(f"  Common ({int(COMMON_THR*100)}-{int(CORE_THR*100)}%): {n_common}")
    print(f"  Intermediate ({int(RARE_THR*100)}-{int(COMMON_THR*100)}%): {n_inter}")
    print(f"  Rare   (< {int(RARE_THR*100)}%) : {n_rare}")

    # ── 5. Presence/absence matrix ─────────────────────────────────────────────
    print("\n=== Building presence/absence matrix ===")
    pa_mat, pa_samples, pa_genes = build_pa_matrix(samples_data)
    write_matrix_csv(pa_mat, pa_samples, pa_genes,
                     os.path.join(OUT_DIR, "presence_absence.csv"))

    # ── 6. Jaccard matrix ──────────────────────────────────────────────────────
    print("\n=== Computing Jaccard similarity matrix ===")
    J = jaccard_matrix(pa_mat)
    write_matrix_csv(np.round(J, 4), pa_samples, pa_samples,
                     os.path.join(OUT_DIR, "jaccard_similarity.csv"))

    # ── 7. Drug class matrix ───────────────────────────────────────────────────
    print("\n=== Building drug class abundance matrix ===")
    dc_mat, dc_samples, dc_classes = drug_class_matrix(samples_data)
    write_matrix_csv(np.round(dc_mat, 2), dc_samples, dc_classes,
                     os.path.join(OUT_DIR, "drug_class_profile.csv"))

    # ── 8. Top lineages ────────────────────────────────────────────────────────
    print("\n=== Computing top ARG-carrying lineages ===")
    lin_rows = top_lineages(samples_data)
    write_csv(lin_rows, os.path.join(OUT_DIR, "top_lineages.csv"))

    # ── 9. Cluster order ───────────────────────────────────────────────────────
    print("\n=== Clustering samples by ARG profile ===")
    order = cluster_order(pa_mat)
    ordered_names = [pa_samples[i] for i in order]
    print(f"  Cluster order (first 5): {ordered_names[:5]}")

    # ── 10. Plots ──────────────────────────────────────────────────────────────
    print("\n=== Generating plots ===")
    plot_gene_prevalence(gf_rows, n_samples,
        os.path.join(OUT_DIR, "01_gene_prevalence.png"))

    plot_drug_class_heatmap(dc_mat, dc_samples, dc_classes, order,
        os.path.join(OUT_DIR, "02_drug_class_heatmap.png"))

    plot_jaccard_heatmap(J, pa_samples, order,
        os.path.join(OUT_DIR, "03_jaccard_heatmap.png"))

    plot_plasmid_burden(stats_rows,
        os.path.join(OUT_DIR, "04_plasmid_burden.png"))

    plot_richness(stats_rows, ordered_names,
        os.path.join(OUT_DIR, "05_sample_richness.png"))

    # ── 11. Console summary ────────────────────────────────────────────────────
    print("\n" + "="*65)
    print("SUMMARY")
    print("="*65)
    print(f"  Samples analysed          : {n_samples}")
    print(f"  Unique ARG genes          : {len(gf_rows)}")

    top5 = gf_rows[:5]
    print(f"\n  Top 5 most prevalent genes:")
    for r in top5:
        print(f"    {r['gene']:20s}  {r['pct_samples']:5.1f}%  [{r['drug_class']}]")

    avg_genes = np.mean([r["n_unique_genes"] for r in stats_rows])
    avg_pf    = np.mean([r["plasmid_frac_copies"] for r in stats_rows])
    print(f"\n  Mean genes per sample     : {avg_genes:.1f}")
    print(f"  Mean plasmid ARG fraction : {avg_pf:.1%}")

    core_genes = [r["gene"] for r in gf_rows if r["prevalence_class"] == "core"]
    if core_genes:
        print(f"\n  Core genes (>= {int(CORE_THR*100)}% of samples):")
        for g in core_genes:
            dc = next(r["drug_class"] for r in gf_rows if r["gene"] == g)
            print(f"    {g:20s}  [{dc}]")
    else:
        print(f"\n  No core genes found (no gene in >= {int(CORE_THR*100)}% of samples).")

    print(f"\n  Output written to: {os.path.abspath(OUT_DIR)}/")
    print("="*65)


if __name__ == "__main__":
    main()
