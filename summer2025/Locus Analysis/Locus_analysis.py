import os
import math
import itertools
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# ====================== Settings ======================
ALPHA = 0.05

# Inputs (exact schemas)
MEAN_CSV      = "mean_test_results/mean_test_summary.csv"                 # Locus,Null_Min,Null_Max,Null_Mean,Test_Statistic,P_Value
CROSS500_CSV  = "crossmatch_results_summary/crossmatch_results_locus_500.csv"    # Locus,nA,nB,a1,Ea1,Va1,dev,pval_exact,pval_normal,elapsed_seconds,timestamp,status,error
CROSS1000_CSV = "crossmatch_results_summary/crossmatch_results_locus_1000.csv"   # same schema

OUTDIR = "locus_analysis_outputs"
os.makedirs(OUTDIR, exist_ok=True)

# Table image settings
TABLE_MAX_ROWS = 60
TABLE_DPI      = 300
# =====================================================

#This script was made to analyze and create graphs for data produced from the crossmatch results for locus trees. Chatgpt was used to help create these functions.

# ------------------------ Loaders ------------------------
def load_mean(path: str) -> pd.DataFrame:
    use_cols = ["Locus","Null_Min","Null_Max","Null_Mean","Test_Statistic","P_Value"]
    df = pd.read_csv(path, usecols=use_cols, encoding="utf-8-sig")
    df = df.dropna(subset=["Locus","P_Value"]).copy()
    df["Locus"]  = df["Locus"].astype(int)
    df["mean_p"] = pd.to_numeric(df["P_Value"], errors="coerce")
    df = df.dropna(subset=["mean_p"])[["Locus","mean_p"]].sort_values("Locus")
    df = df.drop_duplicates(subset=["Locus"], keep="last")
    return df


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

    equal = df["nA"].eq(df["nB"])\

    df["n"] = np.where(equal, df["nA"], np.minimum(df["nA"], df["nB"]))
    df["n"] = df["n"].astype(int)
    if not equal.all():
        print(f"[WARN] {(~equal).sum()} row(s) had nA != nB in {source_label}; used min(nA,nB).")

    df = df[df["n"].isin([50,100,150])].copy()
    df["source"] = source_label
    return df[["Locus","n","pval_exact","pval_normal","source"]].sort_values(["Locus","n"])


# ------------------- Build comparison frames -------------------
def crossmatch_wide(c500: pd.DataFrame, c1000: pd.DataFrame, use_normal=False) -> pd.DataFrame:
    """
    Returns wide table with columns:
      Locus, crossmatch500_n50, crossmatch1000_n50, crossmatch500_n100, crossmatch1000_n100, crossmatch500_n150, crossmatch1000_n150
    """
    val_col = "pval_normal" if use_normal else "pval_exact"
    w500 = c500.pivot(index="Locus", columns="n", values=val_col)
    w1000 = c1000.pivot(index="Locus", columns="n", values=val_col)
    w500 = w500.rename(columns={50:"crossmatch500_n50", 100:"crossmatch500_n100", 150:"crossmatch500_n150"})
    w1000 = w1000.rename(columns={50:"crossmatch1000_n50", 100:"crossmatch1000_n100", 150:"crossmatch1000_n150"})
    wide = (w500.join(w1000, how="inner")).reset_index().sort_values("Locus")
    return wide


# ------------------ Table to PNG ------------------
def save_table_image(df: pd.DataFrame, outpath: str, max_rows=60, dpi=300, rows="head"):
    """
    Render a dataframe to a PNG using matplotlib's table artist.
    - rows="head" -> take the first max_rows after sorting by Locus ascending
    - rows="all"  -> render all rows (may produce a very tall image)
    """
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


# ---------------------- Plot helpers ----------------------
def save_png(outbase):
    """Save only PNG."""
    plt.savefig(f"{outbase}.png", dpi=300, bbox_inches="tight")
    plt.close()


def plot_line_with_mean(sub: pd.DataFrame, n: int, outbase: str,
                        order_by: str = "mean", desc: bool = True):
    """
    order_by: 'mean' | 'cross500' | 'cross1000'
    desc: True = highest p-value first
    """
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


# ---------------------- Binarized agreement toolkit ----------------------
def binarize_tests(df: pd.DataFrame, specs: dict) -> pd.DataFrame:
    """
    Produce a binary matrix of tests (1 = two clusters) based on p-value alpha thresholds (or custom rules).
    """
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


def summarize_binary(bin_df: pd.DataFrame) -> pd.DataFrame:
    cols = list(bin_df.columns)
    pos_counts = bin_df[cols].sum().sort_values(ascending=False)
    totals = pd.Series({c: bin_df[c].notna().sum() for c in cols})
    rate = (pos_counts / totals).rename("PositiveRate")
    return pd.concat([pos_counts.rename("Positives"), totals.rename("N"), rate], axis=1)


# ------------------------- Summary Table Function -------------------------
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


def plot_binary_heatmap(binary_csv, out_png):
    """Create a simple 3-row heatmap from binary_matrix.csv."""
    df = pd.read_csv(binary_csv)
    needed = ["mean_p", "crossmatch500_n50", "crossmatch1000_n50"]
    df = df[["Locus"] + needed].dropna()
    df[needed] = df[needed].astype(int)

    # reshape so rows = test methods, columns = loci
    heatmap_data = df.set_index("Locus")[needed].T


    plt.figure(figsize=(0.25 * len(df) + 3, 2.5))
    sns.heatmap(heatmap_data, cmap="Greys", cbar_kws={"label": "1 = two clusters"})
    plt.title("Binary Cluster Detection (Mean vs Crossmatch)")
    plt.xlabel("Locus")
    plt.ylabel("")
    plt.tight_layout()
    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()


# ------------------------- Main -------------------------
def main():
    # Load data
    mean_df   = load_mean(MEAN_CSV)
    cross500  = load_cross(CROSS500_CSV,  "first/second-500")
    cross1000 = load_cross(CROSS1000_CSV, "all-1000")

    # Build crossmatch wide
    wide = crossmatch_wide(cross500, cross1000, use_normal=False)
    if wide.empty:
        raise SystemExit("No overlapping loci between crossmatch sources.")

    wide = (wide.set_index("Locus")
                .join(mean_df.set_index("Locus")["mean_p"], how="left")
                .reset_index())

    # Summary table (keep PNG + CSV)
    summary_table = create_summary_table(wide)
    summary_csv = os.path.join(OUTDIR, "hypothesis_testing_summary.csv")
    summary_table.to_csv(summary_csv, index=False)

    save_table_image(summary_table, os.path.join(OUTDIR, "hypothesis_testing_summary.png"),
                     max_rows=None, dpi=TABLE_DPI, rows="all")

    print("Summary Table:")
    print(summary_table.to_string(index=False))

    # Save wide CSV (no wide image per request)
    wide_csv = os.path.join(OUTDIR, "wide_with_mean.csv")
    wide.to_csv(wide_csv, index=False)

    # ------------- Phase A: Crossmatch-only, per n -------------
    for n in (50, 100, 150):
        c500_col  = f"crossmatch500_n{n}"
        c1000_col = f"crossmatch1000_n{n}"
        if not all(c in wide.columns for c in [c500_col, c1000_col]):
            continue

        sub = wide[["Locus","mean_p", c500_col, c1000_col]].dropna(subset=[c500_col, c1000_col])

        plot_line_with_mean(sub, n=n, outbase=f"ordered_by_mean_n{n}", order_by="mean")

    # ------------- Phase B: Final three-way @ n=150 -------------
    # Keep CSV 
    need_cols = ["mean_p","crossmatch500_n150","crossmatch1000_n150"]
    if all(c in wide.columns for c in need_cols):
        final150 = wide[["Locus"] + need_cols].dropna().sort_values("Locus")
        final150.to_csv(os.path.join(OUTDIR, "final_compare_n150.csv"), index=False)
    else:
        print("[INFO] Skipping Phase B: missing one of mean_p, crossmatch500_n150, crossmatch1000_n150.")

    # ------------- Binary agreement analysis -------------
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

    # Save binary matrix only (plus a summary CSV for counts/rates)
    bin_df.reset_index().to_csv(os.path.join(OUTDIR, "binary_matrix.csv"), index=False)
    summarize_binary(bin_df).to_csv(os.path.join(OUTDIR, "binary_summary.csv"))

    binary_csv = os.path.join(OUTDIR, "binary_matrix.csv")
    heatmap_png = os.path.join(OUTDIR, "binary_heatmap_3rows.png")
    plot_binary_heatmap(binary_csv, heatmap_png)
    print(f"Saved 3-row binary heatmap to: {heatmap_png}")


if __name__ == "__main__":
    main()
