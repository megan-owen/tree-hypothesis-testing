import os
import math
import itertools
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns

# ====================== Settings ======================
ALPHA = 0.05

# Input file paths
MEAN_CSV      = "batch_output/new_permutation_summary.csv"
CROSS500_CSV  = "crossmatch_results_summary/crossmatch_results_locus_500.csv"
CROSS1000_CSV = "crossmatch_results_summary/crossmatch_results_locus_1000.csv"

# Additional binary test results to merge into heatmap
EXTRA_BIN_CSVS = [
    ("mds_binary", "Locus Analysis/mds.csv"),
    ("odd","Locus Analysis/odd_cases_binary.csv"),
    ("Extreme odd","Locus Analysis/odd_cases_binary.csv"),
]

OUTDIR = "Locus Analysis/locus_analysis_outputs"
os.makedirs(OUTDIR, exist_ok=True)

# Table image settings
TABLE_MAX_ROWS = 60
TABLE_DPI      = 300
# =====================================================

# This script analyzes and visualizes crossmatch test results for locus trees.
# It compares permutation test p-values with crossmatch p-values at different sample sizes.

# ======================== Data Loaders ========================

# Load p-values from permutation test CSV.
# Handles two schemas: (1) filename-based, (2) direct Locus column.
# Auto-detects if header row is missing.
def load_mean(path: str) -> pd.DataFrame:
    # Check if file has proper header
    df_test = pd.read_csv(path, nrows=1, sep=None, engine="python", encoding="utf-8-sig")
    first_col = str(df_test.columns[0]).lower()
    has_header = "filename" in first_col or "locus" in first_col or "rootedtree" not in first_col.lower()
    
    if has_header:
        df = pd.read_csv(path, sep=None, engine="python", encoding="utf-8-sig", on_bad_lines="skip")
    else:
        # No header - manually specify expected columns
        df = pd.read_csv(
            path, 
            header=None,
            names=["filename", "n_perm", "null_min", "null_max", "null_mean", "test_stat", "p_value"],
            sep=None, 
            engine="python", 
            encoding="utf-8-sig", 
            on_bad_lines="skip"
        )

    cols = {c.lower(): c for c in df.columns}
    
    # Schema 1: filename column (e.g., rootedtree_123.txt)
    if "filename" in cols:
        pcol = cols.get("p_value") or cols.get("p value") or cols.get("pvalue")
        if not pcol:
            for c in df.columns:
                if re.fullmatch(r"p[_\s-]*value", c, flags=re.I):
                    pcol = c
                    break
        if not pcol:
            raise ValueError(f"Could not find p_value column. Available: {list(df.columns)}")
        
        tmp = df[[cols["filename"], pcol]].copy()
        tmp.columns = ["filename", "mean_p"]
        # Extract locus number from filename (rootedtree_123.txt -> 123)
        tmp["Locus"] = tmp["filename"].astype(str).str.extract(r"(\d+)").astype(float)
        tmp["mean_p"] = pd.to_numeric(tmp["mean_p"], errors="coerce")
        out = tmp.dropna(subset=["Locus","mean_p"]).copy()
        out["Locus"] = out["Locus"].astype(int)
        return out.sort_values("Locus").drop_duplicates("Locus", keep="last")[["Locus","mean_p"]]
    
    # Schema 2: direct Locus column
    if "locus" in cols:
        pcol = cols.get("p_value") or cols.get("p value") or cols.get("pvalue")
        if pcol is None:
            for c in df.columns:
                if re.fullmatch(r"p[_\s-]*value", c, flags=re.I):
                    pcol = c
                    break
        if pcol:
            out = df[[cols["locus"], pcol]].copy()
            out.columns = ["Locus", "mean_p"]
            out["Locus"]  = pd.to_numeric(out["Locus"], errors="coerce")
            out["mean_p"] = pd.to_numeric(out["mean_p"], errors="coerce")
            out = out.dropna(subset=["Locus","mean_p"])
            out["Locus"] = out["Locus"].astype(int)
            return out.sort_values("Locus").drop_duplicates("Locus", keep="last")[["Locus","mean_p"]]

    raise ValueError(f"Unsupported schema. Expected (Locus, p_value) or (filename, p_value). Found: {list(df.columns)}")


# Load crossmatch summary CSV with p-values per sample size.
# Returns tidy format: Locus, n, pval_exact, pval_normal, source.
def load_cross(path: str, source_label: str) -> pd.DataFrame:
    names = ["Locus","nA","nB","a1","Ea1","Va1","dev",
             "pval_exact","pval_normal","elapsed_seconds","timestamp","status","error"]
    try:
        df = pd.read_csv(path, sep=",", engine="python", encoding="utf-8-sig", skipinitialspace=True)
        if not set(names).issubset(set(df.columns)):
            raise ValueError("Header mismatch")
    except Exception:
        df = pd.read_csv(path, header=None, names=names, sep=",", engine="python",
                         encoding="utf-8-sig", skipinitialspace=True)

    if "status" in df.columns:
        df = df[df["status"].fillna("OK").astype(str).str.upper().eq("OK")]

    for c in ["Locus","nA","nB","pval_exact","pval_normal"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    df = df.dropna(subset=["Locus","nA","nB","pval_exact"]).copy()
    df["Locus"] = df["Locus"].astype(int)

    equal = df["nA"].eq(df["nB"])
    df["n"] = np.where(equal, df["nA"], np.minimum(df["nA"], df["nB"]))
    df["n"] = df["n"].astype(int)
    if not equal.all():
        print(f"[WARN] {(~equal).sum()} row(s) had nA != nB in {source_label}; used min(nA,nB).")

    df = df[df["n"].isin([50,100,150])].copy()
    df["source"] = source_label
    return df[["Locus","n","pval_exact","pval_normal","source"]].sort_values(["Locus","n"])


# Pivot crossmatch results into wide format: Locus, crossmatch500_n50, crossmatch1000_n50, etc.
def crossmatch_wide(c500: pd.DataFrame, c1000: pd.DataFrame, use_normal=False) -> pd.DataFrame:
    val_col = "pval_normal" if use_normal else "pval_exact"
    w500 = c500.pivot(index="Locus", columns="n", values=val_col)
    w1000 = c1000.pivot(index="Locus", columns="n", values=val_col)
    w500 = w500.rename(columns={50:"crossmatch500_n50", 100:"crossmatch500_n100", 150:"crossmatch500_n150"})
    w1000 = w1000.rename(columns={50:"crossmatch1000_n50", 100:"crossmatch1000_n100", 150:"crossmatch1000_n150"})
    wide = (w500.join(w1000, how="inner")).reset_index().sort_values("Locus")
    return wide


# ======================== Visualization Utilities ========================

# Render dataframe as PNG table image.
def save_table_image(df: pd.DataFrame, outpath: str, max_rows=60, dpi=300, rows="head"):
    import matplotlib.pyplot as plt

    data = df.copy()
    if "Locus" in data.columns:
        data = data.sort_values("Locus")

    if rows == "head":
        sub = data if max_rows is None else data.head(max_rows)
    elif rows == "all":
        sub = data
    else:
        sub = data

    nrows, ncols = sub.shape
    cell_w = 1.4
    cell_h = 0.32
    fig_w = max(8, min(28, cell_w * (ncols + 1)))
    fig_h = max(4, min(28, cell_h * (nrows + 2)))

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.axis('off')

    tbl = ax.table(
        cellText=sub.values,
        colLabels=sub.columns,
        cellLoc='center',
        loc='center'
    )
    tbl.auto_set_font_size(False)
    font_size = 9 if ncols <= 8 else 8 if ncols <= 10 else 7
    tbl.set_fontsize(font_size)
    tbl.scale(1, 1.2)

    plt.tight_layout()
    plt.savefig(outpath, dpi=dpi, bbox_inches='tight')
    plt.close()


def save_png(outbase):
    plt.savefig(f"{outbase}.png", dpi=300, bbox_inches="tight")
    plt.close()


# Plot p-values ordered by specified metric with alpha threshold line.
def plot_line_with_mean(sub: pd.DataFrame, n: int, outbase: str,
                        order_by: str = "mean", desc: bool = True):
    aliases = {"orange": "cross500", "green": "cross1000"}
    order_by = aliases.get(order_by, order_by)
    metric_map = {
        "mean": "mean_p",
        "cross500": f"crossmatch500_n{n}",
        "cross1000": f"crossmatch1000_n{n}",
    }
    key = metric_map.get(order_by)
    if key is None or key not in sub.columns:
        raise ValueError(f"Cannot order by '{order_by}': column '{key}' not found.")

    sub = sub.copy()
    sub[key] = pd.to_numeric(sub[key], errors="coerce")
    sub = sub.dropna(subset=[key]).sort_values([key, "Locus"], ascending=[not desc, True])

    x = np.arange(1, len(sub) + 1)
    plt.figure(figsize=(11,5))

    if "mean_p" in sub.columns and sub["mean_p"].notna().any():
        plt.plot(x, sub["mean_p"].values, marker="o", ms=2.5, linestyle="", label="Permutation / Mean test")

    if f"crossmatch500_n{n}" in sub.columns:
        plt.plot(x, sub[f"crossmatch500_n{n}"].values, marker="o", ms=2.5, linestyle="",
                 label=f"Crossmatch (first/second 500), n={n}")

    if f"crossmatch1000_n{n}" in sub.columns:
        plt.plot(x, sub[f"crossmatch1000_n{n}"].values, marker="o", ms=2.5, linestyle="",
                 label=f"Crossmatch (all 1000), n={n}")

    plt.axhline(ALPHA, linestyle="--", label=f"α = {ALPHA}")
    plt.ylabel("p-value")
    plt.xlabel(f"Loci ordered by {order_by} p-value ({'desc' if desc else 'asc'})")

    step = max(1, len(sub)//15)
    plt.xticks(x[::step], sub["Locus"].values[::step])

    plt.title(f"p-values across loci (n={n})")
    plt.legend(ncol=2)
    plt.tight_layout()
    save_png(os.path.join(OUTDIR, outbase))


# ======================== Binary Analysis ========================

# Convert test columns to binary (0/1) based on alpha threshold or custom rules.
def binarize_tests(df: pd.DataFrame, specs: dict) -> pd.DataFrame:
    out = {}
    for col, rule in specs.items():
        if col not in df.columns:
            continue
        s = pd.to_numeric(df[col], errors="coerce")
        t = rule.get("type", "p")
        if t == "p":
            thr = float(rule.get("alpha", ALPHA))
            out[col] = (s < thr).astype(int)
        elif t == "greater":
            thr = float(rule["thr"])
            out[col] = (s > thr).astype(int)
        elif t == "less":
            thr = float(rule["thr"])
            out[col] = (s < thr).astype(int)
        else:
            out[col] = rule["fn"](s).astype(int)
    bin_df = pd.DataFrame(out, index=df.index)
    return bin_df


# Compute summary stats (positives, total, rate) for binary columns.
def summarize_binary(bin_df: pd.DataFrame) -> pd.DataFrame:
    cols = list(bin_df.columns)
    pos_counts = bin_df[cols].sum().sort_values(ascending=False)
    totals = pd.Series({c: bin_df[c].notna().sum() for c in cols})
    rate = (pos_counts / totals).rename("PositiveRate")
    return pd.concat([pos_counts.rename("Positives"), totals.rename("N"), rate], axis=1)


# Generate binary heatmap showing which tests rejected H0 for each locus.
def plot_binary_heatmap(binary_csv, out_png):
    df = pd.read_csv(binary_csv).sort_values("Locus")

    # Exclude n=100 and n=150 from heatmap (focus on n=50)
    EXCLUDE_FOR_HEATMAP = {"crossmatch500_n100", "crossmatch1000_n100","crossmatch500_n150", "crossmatch1000_n150"}
    cols = [c for c in df.columns if c != "Locus" and c not in EXCLUDE_FOR_HEATMAP]

    if not cols:
        print("[WARN] No binary columns found to plot.")
        return

    for c in cols:
        df[c] = pd.to_numeric(df[c], errors="coerce").round().clip(0, 1).fillna(0).astype(int)

    heatmap_data = df.set_index("Locus")[cols].T

    plt.figure(figsize=(0.25 * len(df) + 3, 2 + 0.35 * len(cols)))

    cmap_gray = mcolors.ListedColormap(["white", "0.5"])

    ax = sns.heatmap(
        heatmap_data,
        cmap=cmap_gray,
        cbar=False,
        linewidths=0.5,
        linecolor="black",
        square=True
    )

    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    plt.title("Binary Signals by Locus")
    plt.xlabel("Locus")
    plt.ylabel("")
    plt.tight_layout()
    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()


# Build summary table: rejections, totals, mean p-value per method and sample size.
def create_summary_table(wide: pd.DataFrame) -> pd.DataFrame:
    summary_data = []
    for n in [50, 100, 150]:
        col_500 = f"crossmatch500_n{n}"
        if col_500 in wide.columns:
            pvals_500 = wide[col_500].dropna()
            summary_data.append({
                'Sample_Size': n,
                'Method': 'Crossmatch (500/500)',
                'Rejected': (pvals_500 < ALPHA).sum(),
                'Failed_to_Reject': (pvals_500 >= ALPHA).sum(),
                'Total': len(pvals_500),
                'Mean_p_value': pvals_500.mean()
            })
        col_1000 = f"crossmatch1000_n{n}"
        if col_1000 in wide.columns:
            pvals_1000 = wide[col_1000].dropna()
            summary_data.append({
                'Sample_Size': n,
                'Method': 'Crossmatch (all 1000)',
                'Rejected': (pvals_1000 < ALPHA).sum(),
                'Failed_to_Reject': (pvals_1000 >= ALPHA).sum(),
                'Total': len(pvals_1000),
                'Mean_p_value': pvals_1000.mean()
            })
    if 'mean_p' in wide.columns:
        pvals_mean = wide['mean_p'].dropna()
        summary_data.append({
            'Sample_Size': 'N/A',
            'Method': 'Permutation/Mean Test',
            'Rejected': (pvals_mean < ALPHA).sum(),
            'Failed_to_Reject': (pvals_mean >= ALPHA).sum(),
            'Total': len(pvals_mean),
            'Mean_p_value': pvals_mean.mean()
        })
    return pd.DataFrame(summary_data)


# ======================== Main Pipeline ========================

def main():
    # Load all input data
    mean_df   = load_mean(MEAN_CSV)
    cross500  = load_cross(CROSS500_CSV,  "first/second-500")
    cross1000 = load_cross(CROSS1000_CSV, "all-1000")

    # Merge into wide format
    wide = crossmatch_wide(cross500, cross1000, use_normal=False)
    if wide.empty:
        raise SystemExit("No overlapping loci between crossmatch sources.")

    wide = (wide.set_index("Locus")
                .join(mean_df.set_index("Locus")["mean_p"], how="left")
                .reset_index())

    # Generate and save summary table
    summary_table = create_summary_table(wide)
    summary_csv = os.path.join(OUTDIR, "hypothesis_testing_summary.csv")
    summary_table.to_csv(summary_csv, index=False)
    save_table_image(summary_table, os.path.join(OUTDIR, "hypothesis_testing_summary.png"),
                     max_rows=None, dpi=TABLE_DPI, rows="all")

    print("Summary Table:")
    print(summary_table.to_string(index=False))

    # Save wide CSV
    wide_csv = os.path.join(OUTDIR, "wide_with_mean.csv")
    wide.to_csv(wide_csv, index=False)

    # Create p-value comparison plots for each sample size
    for n in (50, 100, 150):
        c500_col  = f"crossmatch500_n{n}"
        c1000_col = f"crossmatch1000_n{n}"
        if not all(c in wide.columns for c in [c500_col, c1000_col]):
            continue
        sub = wide[["Locus","mean_p", c500_col, c1000_col]].dropna(subset=[c500_col, c1000_col])
        plot_line_with_mean(sub, n=n, outbase=f"ordered_by_mean_n{n}", order_by="mean")

    # Export final comparison at n=150
    need_cols = ["mean_p","crossmatch500_n150","crossmatch1000_n150"]
    if all(c in wide.columns for c in need_cols):
        final150 = wide[["Locus"] + need_cols].dropna().sort_values("Locus")
        final150.to_csv(os.path.join(OUTDIR, "final_compare_n150.csv"), index=False)
    else:
        print("[INFO] Skipping n=150 comparison: missing required columns.")

    # Binary analysis: convert p-values to 0/1 (reject/fail to reject)
    specs = {}
    if "mean_p" in wide.columns:
        specs["mean_p"] = {"type":"p", "alpha": ALPHA}
    for n in (50,100,150):
        for kind in ["500","1000"]:
            col = f"crossmatch{kind}_n{n}"
            if col in wide.columns:
                specs[col] = {"type":"p", "alpha": ALPHA}

    bin_df = binarize_tests(wide.set_index("Locus"), specs)
    bin_df.index.name = "Locus"

    # Merge additional binary test results
    for col_name, path in EXTRA_BIN_CSVS:
        try:
            s = (pd.read_csv(path, encoding="utf-8-sig")
                .rename(columns=lambda c: c.strip()))
            val_col = next(c for c in s.columns if c.strip().lower() != "locus")
            s = s[["Locus", val_col]].dropna()
            s["Locus"] = s["Locus"].astype(int)
            series = pd.to_numeric(s[val_col], errors="coerce").round().clip(0,1).astype("Int64")
            bin_df = bin_df.join(series.rename(col_name), how="left")
        except Exception as e:
            print(f"[WARN] Could not add {col_name} from {path}: {e}")

    # Save binary outputs and generate heatmap
    bin_df.reset_index().to_csv(os.path.join(OUTDIR, "binary_matrix.csv"), index=False)
    summarize_binary(bin_df).to_csv(os.path.join(OUTDIR, "binary_summary.csv"))

    binary_csv = os.path.join(OUTDIR, "binary_matrix.csv")
    heatmap_png = os.path.join(OUTDIR, "binary_heatmap_3rows.png")
    plot_binary_heatmap(binary_csv, heatmap_png)
    print(f"Saved binary heatmap to: {heatmap_png}")


if __name__ == "__main__":
    main()