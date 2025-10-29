# Crossmatch two-sample test over gene trees using a weighted Robinson–Foulds (wRF)
# distance matrix computed via DendroPy.
#
# Workflow:
#   1) Load two sets of Newick trees (FILE_A, FILE_B).
#   2) Build a symmetric pairwise wRF distance matrix on the concatenated trees.
#   3) Run the Crossmatch test and save a readable report.

import dendropy
import numpy as np
import os
from crossmatch_functions import crossmatchtest
from dendropy.calculate import treecompare  # import for RF distances

# =============================================================================
#                         USER-CONFIGURABLE VARIABLES
# =============================================================================

# === Hardcoded Input Files ===
FILE_A = "output/gts_dendropy_CAT_tauAB-10000.0_tauABC-10500.0_tauRoot-11000.0_pAB-10000_pABC-10000_pRoot-10000"
FILE_B = "output/gts_dendropy_CAT_tauAB-10000.0_tauABC-10100.0_tauRoot-11100.0_pAB-10000_pABC-10000_pRoot-10000"
MAX_TREES = 150 # Cap the number of trees taken from each file (None = all lines)
OUTPUT_DIR = "crossmatch_mean_test_results"


# Load trees from a .tre file. 
def load_trees(filename, max_trees=150):
    return dendropy.TreeList.get(path=filename, schema="newick", preserve_underscores=True)[:max_trees]

# Compute a symmetric pairwise weighted RF distance matrix.
# Parameters:
#   trees (dendropy.TreeList): list of trees to compare.
# Returns:
#   np.ndarray: (n x n) matrix where M[i, j] = wRF(t_i, t_j).
def compute_weighted_rf_matrix(trees):
    n = len(trees)
    M = np.zeros((n, n), dtype=float)
    for i in range(n):
        ti = trees[i]
        for j in range(i + 1, n):
            dist = treecompare.weighted_robinson_foulds_distance(ti, trees[j])
            M[i, j] = M[j, i] = float(dist)  # ensure symmetry and float dtype
    return M



# =============================================================================
#                               MAIN RUNNER
# =============================================================================

# Perform the Crossmatch test between two tree samples and write a report.
# Parameters:
#   fileA (str): path to sample A file (one Newick per line).
#   fileB (str): path to sample B file (one Newick per line).
#   out_dir (str): directory to write the report text file.
#   max_trees (int): cap per file when loading trees.
# Returns:
#   dict: {"a1", "Ea1", "Va1", "dev", "pval (exact)", "pval (normal approx)"}.

def run_crossmatch(fileA, fileB, out_dir="mean_crossmatch_results", max_trees=50):
    os.makedirs(out_dir, exist_ok=True)

    # Load and concatenate trees; create labels (0 for A, 1 for B)
    treesA = load_trees(fileA, max_trees)
    treesB = load_trees(fileB, max_trees)
    all_trees = treesA + treesB
    labels = [0] * len(treesA) + [1] * len(treesB)

    print(f"Loaded {len(treesA)} trees from {fileA}")
    print(f"Loaded {len(treesB)} trees from {fileB}")

    print("Computing weighted RF distance matrix...")
    matrix = compute_weighted_rf_matrix(all_trees)

    print("Running crossmatch test...")
    a1, Ea1, Va1, dev, pval_exact, pval_normal = crossmatchtest(labels, matrix)

    results = {
        "a1": a1,
        "Ea1": Ea1,
        "Va1": Va1,
        "dev": dev,
        "pval (exact)": pval_exact,
        "pval (normal approx)": pval_normal,
    }

    # Construct output filename from the input names
    base1 = os.path.basename(fileA).replace(".tre", "")
    base2 = os.path.basename(fileB).replace(".tre", "")
    out_file = os.path.join(out_dir, f"crossmatch_{base1}_vs_{base2}.txt")

    print(f"\nResults saved to: {out_file}")
    with open(out_file, "w") as f:
        for k, v in results.items():
            print(f"{k}: {v}")
            f.write(f"{k}: {v}\n")

    return results



# === Run test ===
if __name__ == "__main__":
    run_crossmatch(FILE_A, FILE_B, out_dir=OUTPUT_DIR, max_trees=MAX_TREES)