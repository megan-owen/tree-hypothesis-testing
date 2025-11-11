# A general, single-script workflow to run the Crossmatch test using either
#  - weighted Robinson–Foulds (wRF) distances (via DendroPy), or
#  - BHV distances (via gtp.jar + JPype).
#
# All configuration lives at the top; you can also override via CLI (argparse).

import os
import sys
import csv
import time
import random
import argparse
from datetime import datetime

import numpy as np
import dendropy
from dendropy.calculate import treecompare
import jpype
import jpype.imports

from crossmatch_functions import crossmatchtest

# =============================================================================
#                         USER-CONFIGURABLE VARIABLES
# =============================================================================

DEBUG = True                       # If True, print sanity checks and extra logging
METHOD = "wrf"                     # "wrf" or "bhv"
FILE_A = "path/to/sampleA.newick"  # File containing one Newick tree per line
FILE_B = "path/to/sampleB.newick"  # File containing one Newick tree per line

SAMPLE_SIZE_A = 150                # Number of trees sampled from FILE_A
SAMPLE_SIZE_B = 150                # Number of trees sampled from FILE_B
MAX_TREES = None                   # Optional cap per file (None = all)
ROOTED = True                      # Treat trees as rooted (affects BHV)
GTP_JAR_PATH = "gtp.jar"           # Path to gtp.jar (needed only for METHOD="bhv")

OUTPUT_DIR = "crossmatch_results"  # Directory to save outputs
WRITE_CSV_SUMMARY = True           # Write a CSV with one-row summary
SUMMARY_CSV_PATH = os.path.join(OUTPUT_DIR, "crossmatch_summary.csv")

# RNG seeding: if None, seed from time (ms). Set to an int for reproducibility.
SEED = None


# -----------------------------------------------------------------------------
# What you can change (quick guide):
#   - DEBUG:     toggle extra prints and sanity checks
#   - METHOD:    "wrf" (DendroPy weighted RF) or "bhv" (BHV via gtp.jar)
#   - FILE_A/B:  any two files with one Newick tree per line
#   - SAMPLE_SIZE_A/B: per-group sample sizes
#   - MAX_TREES: limit how many lines to read from each file (None = all)
#   - ROOTED:    BHV rootedness flag
#   - GTP_JAR_PATH: where gtp.jar lives (BHV only)
#   - OUTPUT_DIR / SUMMARY_CSV_PATH / WRITE_CSV_SUMMARY: output behavior
# =============================================================================

# =========================
# File & tree utilities
# =========================

# Read a file of Newick strings (one per line).
# Parameters:
#   path (str): path to a text file with one Newick tree per line.
#   max_trees (int|None): optional limit of how many lines to read; None = all.
# Returns:
#   list[str]: list of Newick strings (non-empty, stripped).
def read_newick_file(path, max_trees=None):
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Missing file: {path}")
    with open(path, "r") as f:
        trees = [line.strip() for line in f if line.strip()]
    return trees[:max_trees] if max_trees else trees


# Convert list of Newick strings into a DendroPy TreeList.
# Parameters:
#   lines (list[str]): Newick strings.
# Returns:
#   dendropy.TreeList: list of parsed trees.
def load_treelist_from_lines(lines):
    if dendropy is None:
        raise ImportError("DendroPy is not available; cannot compute wRF distances.")
    tl = dendropy.TreeList()
    for t_str in lines:
        tl.append(dendropy.Tree.get(data=t_str, schema="newick", preserve_underscores=True))
    return tl


# =========================
# Distance computations
# =========================

# Compute weighted Robinson–Foulds distance matrix using DendroPy.
# Parameters:
#   trees (dendropy.TreeList): trees to compare.
# Returns:
#   np.ndarray: symmetric matrix M (n x n) with wRF distances.
def compute_wrf_distance_matrix(trees):
    n = len(trees)
    M = np.zeros((n, n), dtype=float)
    for i in range(n):
        ti = trees[i]
        for j in range(i + 1, n):
            dist = treecompare.weighted_robinson_foulds_distance(ti, trees[j])
            M[i, j] = M[j, i] = dist
        if DEBUG and (i + 1) % 25 == 0:
            print(f"[wrf] row {i+1}/{n} done")
    return M


# Start the JVM for BHV calculations if not already started, and bind classes.
# Parameters:
#   jar_path (str): path to gtp.jar.
# Returns:
#   tuple(PolyMain, PhyloTree): Java classes for geodesic and tree objects.
def ensure_jvm_and_bind_classes(jar_path):
    if not os.path.isfile(jar_path):
        raise FileNotFoundError(f"gtp.jar not found at: {jar_path}")
    if not jpype.isJVMStarted():
        jpype.startJVM(classpath=[jar_path])
    PolyMain = jpype.JClass("polyAlg.PolyMain")
    PhyloTree = jpype.JClass("distanceAlg1.PhyloTree")
    return PolyMain, PhyloTree


# Compute BHV pairwise distance matrix using Java classes.
# Parameters:
#   newick_strings (list[str]): concatenated [A + B] Newick strings.
#   rooted (bool): whether to treat trees as rooted in the BHV space.
#   classes (tuple): (PolyMain, PhyloTree) returned by ensure_jvm_and_bind_classes.
# Returns:
#   np.ndarray: symmetric matrix M (n x n) with BHV distances (float).
def compute_bhv_distance_matrix(newick_strings, rooted, classes):
    PolyMain, PhyloTree = classes
    n = len(newick_strings)
    M = np.zeros((n, n), dtype=float)

    print("[bhv] Computing pairwise BHV distances...")
    for i in range(n):
        for j in range(i + 1, n):
            try:
                t1 = PhyloTree(newick_strings[i], rooted)
                t2 = PhyloTree(newick_strings[j], rooted)
                dist = PolyMain.getGeodesic(t1, t2, None).getDist()
                M[i, j] = M[j, i] = float(dist)
            except Exception as e:
                if DEBUG:
                    print(f"[bhv] Error between {i} and {j}: {e}")
                M[i, j] = M[j, i] = np.nan
        if DEBUG and (i + 1) % 10 == 0:
            print(f"[bhv] row {i+1}/{n} done")
    return M


# =========================
# Sampling & sanity checks
# =========================

# Select two disjoint random index sets (A and B) from 0..(N-1) independently per file.
# Parameters:
#   nA_avail (int): total available in file A.
#   nB_avail (int): total available in file B.
#   nA (int): number to sample from A.
#   nB (int): number to sample from B.
# Returns:
#   tuple(list[int], list[int]): indicesA, indicesB
def select_indices_two_files(nA_avail, nB_avail, nA, nB):
    if nA > nA_avail:
        raise ValueError(f"Requested nA={nA} exceeds available A={nA_avail}")
    if nB > nB_avail:
        raise ValueError(f"Requested nB={nB} exceeds available B={nB_avail}")
    idxA = random.sample(range(nA_avail), nA)
    idxB = random.sample(range(nB_avail), nB)
    return idxA, idxB


# Basic sanity checks on the distance matrix and labels.
# Parameters:
#   M (np.ndarray): distance matrix (n x n).
#   labels (list[int]): group labels (0 for A, 1 for B), length n.
# Returns:
#   None (raises AssertionError in DEBUG mode).
def sanity_check_matrix_and_labels(M, labels):
    if not DEBUG:
        return
    n = len(labels)
    assert M.shape == (n, n), f"Matrix shape {M.shape} != ({n},{n})"
    assert np.allclose(M, M.T, equal_nan=True), "Matrix not symmetric."
    if np.isnan(M).any():
        print(f"[warn] Distance matrix contains NaNs; crossmatchtest may fail.")


# =========================
# I/O helpers
# =========================

# Append a one-row summary to CSV (header auto-written if missing).
# Parameters:
#   row (dict): keys become columns.
#   csv_path (str): where to write/append.
# Returns:
#   None (writes to file).
def append_summary_row(row, csv_path):
    os.makedirs(os.path.dirname(csv_path), exist_ok=True)
    file_exists = os.path.isfile(csv_path)
    with open(csv_path, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(row.keys()))
        if not file_exists:
            writer.writeheader()
        writer.writerow(row)


# =========================
# Main crossmatch runner
# =========================

# Run a general crossmatch experiment with the chosen method (wRF or BHV).
# Parameters:
#   fileA (str): path to sample A trees (one Newick per line).
#   fileB (str): path to sample B trees (one Newick per line).
#   method (str): "wrf" or "bhv".
#   nA (int): sample size for group A.
#   nB (int): sample size for group B.
#   rooted (bool): used by BHV construction (ignored for wRF).
#   max_trees (int|None): optional limit per file; None = all.
#   out_dir (str): output directory for text/CSV artifacts.
#   gt_jar (str): path to gtp.jar (BHV only).
# Returns:
#   dict: results including stats and p-values.
def run_crossmatch_general(fileA, fileB, method, nA, nB, rooted, max_trees, out_dir, gt_jar):
    os.makedirs(out_dir, exist_ok=True)

    # Load trees
    A_lines = read_newick_file(fileA, max_trees=max_trees)
    B_lines = read_newick_file(fileB, max_trees=max_trees)
    if DEBUG:
        print(f"[debug] Loaded {len(A_lines)} trees from A: {fileA}")
        print(f"[debug] Loaded {len(B_lines)} trees from B: {fileB}")

    # Sample indices independently from each file
    idxA, idxB = select_indices_two_files(len(A_lines), len(B_lines), nA, nB)
    treesA_newick = [A_lines[i] for i in idxA]
    treesB_newick = [B_lines[i] for i in idxB]

    labels = [0] * len(treesA_newick) + [1] * len(treesB_newick)

    # Compute distances
    t0 = time.time()
    if method.lower() == "wrf":
        treesA = load_treelist_from_lines(treesA_newick)
        treesB = load_treelist_from_lines(treesB_newick)
        all_trees = dendropy.TreeList(treesA + treesB)
        M = compute_wrf_distance_matrix(all_trees)
    elif method.lower() == "bhv":
        classes = ensure_jvm_and_bind_classes(gt_jar)
        all_newick = treesA_newick + treesB_newick
        M = compute_bhv_distance_matrix(all_newick, rooted, classes)
    else:
        raise ValueError(f"Unknown method: {method}. Use 'wrf' or 'bhv'.")
    elapsed_dist = time.time() - t0

    # Sanity + test
    sanity_check_matrix_and_labels(M, labels)
    if DEBUG:
        print("[debug] Running crossmatchtest...")
    a1, Ea1, Va1, dev, pval_exact, pval_normal = crossmatchtest(labels, M)

    # Prepare outputs
    results = {
        "fileA": os.path.basename(fileA),
        "fileB": os.path.basename(fileB),
        "method": method.lower(),
        "nA": len(treesA_newick),
        "nB": len(treesB_newick),
        "rooted": rooted if method.lower() == "bhv" else None,
        "a1": a1,
        "Ea1": Ea1,
        "Va1": Va1,
        "dev": dev,
        "pval_exact": pval_exact,
        "pval_normal": pval_normal,
        "elapsed_distance_sec": round(elapsed_dist, 3),
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "status": "OK",
        "error": "",
    }

    # Write a per-run text file
    tag = f"{method}_A{results['nA']}_B{results['nB']}_{int(time.time())}"
    out_txt = os.path.join(out_dir, f"crossmatch_{tag}.txt")
    with open(out_txt, "w") as f:
        for k, v in results.items():
            f.write(f"{k}: {v}\n")

    print("\n=== Crossmatch Results ===")
    for k in ("method", "nA", "nB", "a1", "Ea1", "Va1", "dev", "pval_exact", "pval_normal", "elapsed_distance_sec"):
        print(f"{k}: {results[k]}")
    print(f"Saved: {out_txt}")

    # Optional CSV summary
    if WRITE_CSV_SUMMARY:
        append_summary_row(results, SUMMARY_CSV_PATH)
        if DEBUG:
            print(f"[debug] Appended summary to {SUMMARY_CSV_PATH}")

    return results


# =========================
# Main Execution
# =========================

if __name__ == "__main__":
    # Global RNG seed from current time in milliseconds (different each run)
    random.seed(int(time.time() * 1000))
    print(f"[info] METHOD={METHOD} nA={SAMPLE_SIZE_A} nB={SAMPLE_SIZE_B} ROOTED={ROOTED} DEBUG={DEBUG}")

    t_start = time.time()
    try:
        run_crossmatch_general(
            fileA=FILE_A,
            fileB=FILE_B,
            method=METHOD,
            nA=SAMPLE_SIZE_A,
            nB=SAMPLE_SIZE_B,
            rooted=ROOTED,
            max_trees=MAX_TREES,
            out_dir=OUTPUT_DIR,
            gt_jar=GTP_JAR_PATH,
        )
    finally:
        # Cleanly shutdown JVM if we started it (BHV case)
        if jpype.isJVMStarted():
            try:
                jpype.shutdownJVM()
            except Exception:
                pass
        print(f"\nTotal time: {time.time() - t_start:.2f} s")