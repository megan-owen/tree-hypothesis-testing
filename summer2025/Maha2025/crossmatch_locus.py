import os
import time
import csv
from datetime import datetime
import random

import numpy as np
import dendropy
from dendropy.calculate import treecompare
from crossmatch_functions import crossmatchtest

# =========================
# Configuration & Paths
# =========================

ROOTED_DIR = "rootedtrees"                 # Directory with rooted gene tree files
OUTPUT_DIR = "crossmatch_results_locus"    # Per-locus crossmatch result files

SAVE_DIR = "crossmatch_results_summary"    # Directory for summary CSV
os.makedirs(SAVE_DIR, exist_ok=True)
OUTPUT_CSV = os.path.join(SAVE_DIR, "crossmatch_summary.csv")

BATCH_FSYNC = 10                           # Flush data every 10 loci to avoid data loss

# Split strategy for forming groups from a single locus file:
#   "whole"  -> sample both groups randomly (disjoint) from the entire set
#   "halves" -> split the set into two equal halves (by index); sample A from first half, B from second half
SPLIT_MODE = "halves"                      # "whole" or "halves"

# When sampling from each group
SAMPLE_SIZE_A = 150
SAMPLE_SIZE_B = 150

# Seed RNG using current time in milliseconds
random.seed(int(time.time() * 1000))

# Columns for the summary CSV
CSV_COLUMNS = [
    "Locus",
    "nA", "nB",
    "a1", "Ea1", "Va1", "dev",
    "pval_exact", "pval_normal",
    "elapsed_seconds",
    "timestamp", "status", "error",
]


# =========================
# Helper Functions
# =========================

# Append one row of results to the summary CSV (writes header if missing).
# Parameters:
#   row_dict (dict): values for a single locus run
#   csv_path (str): path to the CSV file
# Returns:
#   None (writes to disk)
def _append_csv_row(row_dict, csv_path=OUTPUT_CSV):
    file_exists = os.path.isfile(csv_path)
    with open(csv_path, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
        if not file_exists:
            writer.writeheader()
        writer.writerow(row_dict)
    # fsync handled in the main loop to reduce overhead


# Convert list of Newick strings into a DendroPy TreeList.
# Parameters:
#   lines (list[str]): each string is a Newick tree
# Returns:
#   dendropy.TreeList
def load_trees_from_lines(lines):
    tl = dendropy.TreeList()
    for t_str in lines:
        tl.append(dendropy.Tree.get(data=t_str, schema="newick", preserve_underscores=True))
    return tl


# Compute weighted Robinson–Foulds distance matrix between all pairs of trees.
# Parameters:
#   trees (dendropy.TreeList)
# Returns:
#   np.ndarray (n x n) symmetric matrix of distances
def compute_weighted_rf_matrix(trees):
    n = len(trees)
    matrix = np.zeros((n, n), dtype=float)
    for i in range(n):
        ti = trees[i]
        for j in range(i + 1, n):
            dist = treecompare.weighted_robinson_foulds_distance(ti, trees[j])
            matrix[i, j] = dist
            matrix[j, i] = dist
    return matrix



# Choose indices for group A and B based on SPLIT_MODE.
# Parameters:
#   tree_lines (list[str]): all trees read from the file
#   nA (int): sample size for group A
#   nB (int): sample size for group B
# Returns:
#   (list[int], list[int]): indicesA, indicesB (disjoint)
def select_indices_for_locus(tree_lines, nA, nB):
    N = len(tree_lines)
    if N < nA + nB and SPLIT_MODE.lower() == "whole":
        raise ValueError(f"Not enough trees for whole-sample draw: need {nA+nB}, have {N}")

    if SPLIT_MODE.lower() == "halves":
        # split point is dynamic based on file length
        mid = N // 2
        first_half = list(range(0, mid))
        second_half = list(range(mid, N))
        if nA > len(first_half):
            raise ValueError(f"nA={nA} exceeds first-half size={len(first_half)} (N={N})")
        if nB > len(second_half):
            raise ValueError(f"nB={nB} exceeds second-half size={len(second_half)} (N={N})")
        indicesA = random.sample(first_half, nA)
        indicesB = random.sample(second_half, nB)

    elif SPLIT_MODE.lower() == "whole":
        all_indices = list(range(N))
        indicesA = random.sample(all_indices, nA)
        remaining = [i for i in all_indices if i not in indicesA]
        if len(remaining) < nB:
            raise ValueError(f"Not enough remaining trees for group B: need {nB}, have {len(remaining)}")
        indicesB = random.sample(remaining, nB)

    else:
        raise ValueError(f"Unknown SPLIT_MODE: {SPLIT_MODE} (use 'whole' or 'halves')")

    # sanity: disjoint
    overlap = set(indicesA) & set(indicesB)
    print(f"[split:{SPLIT_MODE}] N={N} mid={N//2} nA={len(indicesA)} nB={len(indicesB)} overlap={len(overlap)}")
    assert not overlap, "Split overlap detected!"
    return indicesA, indicesB


# Run the Crossmatch test for one locus and write a small .txt file.
# Parameters:
#   treesA (dendropy.TreeList), treesB (dendropy.TreeList)
#   out_dir (str): where to write per-locus txt
#   locus_id (int)
# Returns:
#   dict with stats and p-values
def run_crossmatch(treesA, treesB, out_dir, locus_id):
    os.makedirs(out_dir, exist_ok=True)

    all_trees = treesA + treesB
    assert len(all_trees) == len(treesA) + len(treesB)
    labels = [0] * len(treesA) + [1] * len(treesB)

    print(f"  Computing weighted RF matrix for Locus {locus_id} "
          f"({len(treesA)} A + {len(treesB)} B)...")
    matrix = compute_weighted_rf_matrix(all_trees)

    print(f"  Running crossmatch test for Locus {locus_id}...")
    a1, Ea1, Va1, dev, pval_exact, pval_normal = crossmatchtest(labels, matrix)

    results = {
        "a1": a1,
        "Ea1": Ea1,
        "Va1": Va1,
        "dev": dev,
        "pval_exact": pval_exact,
        "pval_normal": pval_normal,
    }

    out_file = os.path.join(out_dir, f"crossmatch_locus_{locus_id}.txt")
    with open(out_file, "w") as f:
        for k, v in results.items():
            f.write(f"{k}: {v}\n")

    return results


# =========================
# Main Execution
# =========================
if __name__ == "__main__":
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    loci_processed_since_fsync = 0

    # Loop over loci
    for locus_id in range(19, 21):
        file_path = os.path.join(ROOTED_DIR, f"rootedtree_{locus_id}.txt")
        start_wall = time.time()

        # Prepare default row for summary CSV
        row = {
            "Locus": locus_id,
            "nA": 0, "nB": 0,
            "a1": "", "Ea1": "", "Va1": "", "dev": "",
            "pval_exact": "", "pval_normal": "",
            "elapsed_seconds": "",
            "timestamp": datetime.now().isoformat(timespec="seconds"),
            "status": "OK", "error": "",
        }

        print(f"\n=== Locus {locus_id}: {os.path.basename(file_path)} ===")

        try:
            # Load trees for this locus
            if not os.path.isfile(file_path):
                raise FileNotFoundError(f"Missing file: {file_path}")
            with open(file_path, "r") as f:
                tree_lines = [line.strip() for line in f if line.strip()]

            # Choose split strategy ("whole" vs "halves") by SPLIT_MODE
            indicesA, indicesB = select_indices_for_locus(
                tree_lines,
                nA=SAMPLE_SIZE_A,
                nB=SAMPLE_SIZE_B
            )

            # Build TreeLists
            treesA = load_trees_from_lines([tree_lines[i] for i in indicesA])
            treesB = load_trees_from_lines([tree_lines[i] for i in indicesB])

            row["nA"] = len(treesA)
            row["nB"] = len(treesB)

            # Run Crossmatch
            results = run_crossmatch(treesA, treesB, out_dir=OUTPUT_DIR, locus_id=locus_id)

            # Merge results
            row.update({
                "a1": results["a1"],
                "Ea1": results["Ea1"],
                "Va1": results["Va1"],
                "dev": results["dev"],
                "pval_exact": results["pval_exact"],
                "pval_normal": results["pval_normal"],
            })

        except Exception as e:
            row["status"] = "ERROR"
            row["error"] = str(e)
            print(f"  [ERROR] Locus {locus_id}: {e}")

        finally:
            # Log to CSV
            row["elapsed_seconds"] = round(time.time() - start_wall, 3)
            _append_csv_row(row, OUTPUT_CSV)
            loci_processed_since_fsync += 1

            # Periodically flush results
            if loci_processed_since_fsync % BATCH_FSYNC == 0 or locus_id in (152, 304):
                with open(OUTPUT_CSV, "a") as f:
                    f.flush()
                    os.fsync(f.fileno())
                print(f"  [fsync] Flushed CSV after locus {locus_id}")

            # Progress
            if row["status"] == "OK":
                print(f"  Done Locus {locus_id}: "
                      f"a1={row['a1']}, p_exact={row['pval_exact']}, "
                      f"t={row['elapsed_seconds']}s")
            else:
                print(f"  Logged ERROR for Locus {locus_id} (t={row['elapsed_seconds']}s)")

    print(f"\n=== Finished. Summary CSV at: {OUTPUT_CSV} ===")