import os
import time
import csv
from datetime import datetime
import random

import numpy as np
import dendropy
from dendropy.calculate import treecompare
from crossmatch_functions import crossmatchtest

# Configuration & Paths

ROOTED_DIR = "rootedtrees"                # Directory with rooted gene tree files
OUTPUT_DIR = "crossmatch_results_locus"   # Per-locus crossmatch result files

SAVE_DIR = "crossmatch_results_summary"   # Directory for summary CSV
os.makedirs(SAVE_DIR, exist_ok=True)
OUTPUT_CSV = os.path.join(SAVE_DIR, "crossmatch_summary.csv")

BATCH_FSYNC = 10  # Flush data every 10 loci to avoid data loss

# Columns for the summary CSV
CSV_COLUMNS = [
    "Locus",
    "nA", "nB",
    "a1", "Ea1", "Va1", "dev",
    "pval_exact", "pval_normal",
    "elapsed_seconds",
    "timestamp", "status", "error",
]


# Helper Functions

def _append_csv_row(row_dict, csv_path=OUTPUT_CSV):
    """Appends one row of results to the summary CSV (writes header if missing)."""
    file_exists = os.path.isfile(csv_path)
    with open(csv_path, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
        if not file_exists:
            writer.writeheader()
        writer.writerow(row_dict)
        # fsync handled later to reduce overhead.


def load_trees_from_lines(lines):
    """Convert list of Newick strings into a DendroPy TreeList."""
    tl = dendropy.TreeList()
    for t_str in lines:
        tl.append(dendropy.Tree.get(data=t_str, schema="newick", preserve_underscores=True))
    return tl


def compute_weighted_rf_matrix(trees):
    """Compute weighted Robinson-Foulds distance matrix between all pairs of trees."""
    n = len(trees)
    matrix = np.zeros((n, n), dtype=float)
    for i in range(n):
        ti = trees[i]
        for j in range(i + 1, n):
            dist = treecompare.weighted_robinson_foulds_distance(ti, trees[j])
            matrix[i, j] = dist
            matrix[j, i] = dist
    return matrix


def run_crossmatch(treesA, treesB, out_dir, locus_id):
    """Run Crossmatch test for one locus and save results to a .txt file."""
    os.makedirs(out_dir, exist_ok=True)

    # Combine both groups
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

    # Save results to text file
    out_file = os.path.join(out_dir, f"crossmatch_locus_{locus_id}.txt")
    with open(out_file, "w") as f:
        for k, v in results.items():
            f.write(f"{k}: {v}\n")

    return results


# Main Execution
if __name__ == "__main__":
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    random.seed(42)  # Reproducible random splits

    loci_processed_since_fsync = 0

    # Loop over loci
    for locus_id in range(19, 21):
        file_path = os.path.join(ROOTED_DIR, f"rootedtree_{locus_id}.txt")
        start_wall = time.time()

        # Initialize default row for summary CSV
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
            # Step 1: Load all trees for this locus
            if not os.path.isfile(file_path):
                raise FileNotFoundError(f"Missing file: {file_path}")

            with open(file_path, "r") as f:
                tree_lines = [line.strip() for line in f if line.strip()]

            # OPTION 1 — Fixed halves (first 500 vs second 500 trees)
            '''
            # This is the standard configuration for comparing.

            if len(tree_lines) < 1000:
                raise ValueError(f"Expected 1000 trees, got {len(tree_lines)}")

            first_half = list(range(0, 500))
            second_half = list(range(500, 1000))

            # Randomly select 150 from each half
            indicesA = random.sample(first_half, 150)
            indicesB = random.sample(second_half, 150)
            '''

            # OPTION 2 — Random split across all 1000 trees 
            all_indices = list(range(len(tree_lines)))
            indicesA = random.sample(all_indices, 150)

            remaining_indices = [i for i in all_indices if i not in indicesA]
            if len(remaining_indices) < 151:
                raise ValueError(f"Not enough remaining trees for group B: {len(remaining_indices)}")
            indicesB = random.sample(remaining_indices, 150)

            # Sanity check: ensure disjoint groups
            overlap = set(indicesA) & set(indicesB)
            print(f"[check] locus {locus_id}: nA={len(indicesA)} nB={len(indicesB)} overlap={len(overlap)}")
            assert not overlap, "Random split overlap detected!"

            # Step 2: Load trees into DendroPy TreeLists
            treesA = load_trees_from_lines([tree_lines[i] for i in indicesA])
            treesB = load_trees_from_lines([tree_lines[i] for i in indicesB])

            row["nA"] = len(treesA)
            row["nB"] = len(treesB)

            # Step 3: Run Crossmatch test

            results = run_crossmatch(treesA, treesB, out_dir=OUTPUT_DIR, locus_id=locus_id)

            # Merge results into row
            row.update({
                "a1": results["a1"],
                "Ea1": results["Ea1"],
                "Va1": results["Va1"],
                "dev": results["dev"],
                "pval_exact": results["pval_exact"],
                "pval_normal": results["pval_normal"],
            })

        except Exception as e:
            # Record any errors (missing file, invalid data, etc.)
            row["status"] = "ERROR"
            row["error"] = str(e)
            print(f"  [ERROR] Locus {locus_id}: {e}")

        finally:

            # Step 4: Log results to summary CSV
            
            row["elapsed_seconds"] = round(time.time() - start_wall, 3)
            _append_csv_row(row, OUTPUT_CSV)
            loci_processed_since_fsync += 1

            # Periodically flush results to disk
            if loci_processed_since_fsync % BATCH_FSYNC == 0 or locus_id in (152, 304):
                with open(OUTPUT_CSV, "a") as f:
                    f.flush()
                    os.fsync(f.fileno())
                print(f"  [fsync] Flushed CSV after locus {locus_id}")

            # Progress log
            if row["status"] == "OK":
                print(f"  Done Locus {locus_id}: "
                      f"a1={row['a1']}, p_exact={row['pval_exact']}, "
                      f"t={row['elapsed_seconds']}s")
            else:
                print(f"  Logged ERROR for Locus {locus_id} (t={row['elapsed_seconds']}s)")

    print(f"\n=== Finished. Summary CSV at: {OUTPUT_CSV} ===")
