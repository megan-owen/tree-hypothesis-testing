import os
import time
import csv
from datetime import datetime

import numpy as np
import dendropy
from dendropy.calculate import treecompare
from crossmatch_functions import crossmatchtest

ROOTED_DIR = "rootedtrees"
OUTPUT_DIR = "crossmatch_results_locus"

# CSV settings
SAVE_DIR = "crossmatch_results_summary"
os.makedirs(SAVE_DIR, exist_ok=True)
OUTPUT_CSV = os.path.join(SAVE_DIR, "crossmatch_summary.csv")
BATCH_FSYNC = 10  # force flush/fsync every 10 loci

CSV_COLUMNS = [
    "Locus",
    "nA",
    "nB",
    "a1",
    "Ea1",
    "Va1",
    "dev",
    "pval_exact",
    "pval_normal",
    "elapsed_seconds",
    "timestamp",
    "status",
    "error",
]

def _append_csv_row(row_dict, csv_path=OUTPUT_CSV):
    file_exists = os.path.isfile(csv_path)
    with open(csv_path, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
        if not file_exists:
            writer.writeheader()
        writer.writerow(row_dict)
        # caller is responsible for batching fsyncs; this keeps per-row overhead low

def load_trees_from_lines(lines):
    """Convert a list of Newick strings into a DendroPy TreeList."""
    tl = dendropy.TreeList()
    for t_str in lines:
        tl.append(dendropy.Tree.get(data=t_str, schema="newick", preserve_underscores=True))
    return tl

def compute_weighted_rf_matrix(trees):
    """Compute symmetric weighted RF matrix for a list/TreeList of trees."""
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
    os.makedirs(out_dir, exist_ok=True)

    all_trees = treesA + treesB
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

if __name__ == "__main__":
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    loci_processed_since_fsync = 0
    for locus_id in range(152, 305):  
        file_path = os.path.join(ROOTED_DIR, f"rootedtree_{locus_id}.txt")
        start_wall = time.time()
        row = {
            "Locus": locus_id,
            "nA": 0,
            "nB": 0,
            "a1": "",
            "Ea1": "",
            "Va1": "",
            "dev": "",
            "pval_exact": "",
            "pval_normal": "",
            "elapsed_seconds": "",
            "timestamp": datetime.now().isoformat(timespec="seconds"),
            "status": "OK",
            "error": "",
        }

        print(f"\n=== Locus {locus_id}: {os.path.basename(file_path)} ===")
        try:
            if not os.path.isfile(file_path):
                raise FileNotFoundError(f"Missing file: {file_path}")

            with open(file_path, "r") as f:
                tree_lines = [line.strip() for line in f if line.strip()]

            if len(tree_lines) < 551:
                raise ValueError(f"Expected at least 551 lines; got {len(tree_lines)}")

            treesA_lines = tree_lines[:50]
            treesB_lines = tree_lines[500:550]

            treesA = load_trees_from_lines(treesA_lines)
            treesB = load_trees_from_lines(treesB_lines)

            row["nA"] = len(treesA)
            row["nB"] = len(treesB)

            results = run_crossmatch(treesA, treesB, out_dir=OUTPUT_DIR, locus_id=locus_id)

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
            row["elapsed_seconds"] = round(time.time() - start_wall, 3)
            _append_csv_row(row, OUTPUT_CSV)
            loci_processed_since_fsync += 1

            # Force a flush & fsync every BATCH_FSYNC loci
            if loci_processed_since_fsync % BATCH_FSYNC == 0 or locus_id in (152, 304):
                with open(OUTPUT_CSV, "a") as f:
                    f.flush()
                    os.fsync(f.fileno())
                print(f"  [fsync] Flushed CSV after locus {locus_id}")

            # Small progress line
            if row["status"] == "OK":
                print(f"  Done Locus {locus_id}: "
                      f"a1={row['a1']}, p_exact={row['pval_exact']}, "
                      f"t={row['elapsed_seconds']}s")
            else:
                print(f"  Logged ERROR for Locus {locus_id} (t={row['elapsed_seconds']}s)")

    print(f"\n=== Finished. Summary CSV at: {OUTPUT_CSV} ===")
