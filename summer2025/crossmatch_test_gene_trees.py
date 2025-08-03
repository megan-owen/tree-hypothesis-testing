import dendropy
import numpy as np
import os
from crossmatch_functions import crossmatchtest
from dendropy.calculate import treecompare  # import for RF distances


# === Hardcoded Input Files ===
FILE_A = "output/gts_dendropy_CAT_tauAB-10000.0_tauABC-10500.0_tauRoot-11000.0_pAB-10000_pABC-10000_pRoot-10000"
FILE_B = "output/gts_dendropy_CAT_tauAB-10000.0_tauABC-10100.0_tauRoot-11100.0_pAB-10000_pABC-10000_pRoot-10000"
MAX_TREES = 165
OUTPUT_DIR = "crossmatch_results"


# Load trees from a .tre file. Temporarily testing with only 50 trees from each sample.
def load_trees(filename, max_trees=165):
    return dendropy.TreeList.get(path=filename, schema="newick", preserve_underscores=True)[:max_trees]

# Compute pairwise weighted RF distance matrix
# From deepseek
def compute_weighted_rf_matrix(trees):
    n = len(trees)
    matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            dist = treecompare.weighted_robinson_foulds_distance(trees[i], trees[j])
            matrix[i][j] = matrix[j][i] = dist  # symmetric matrix
    return matrix

# Perform the cross-match test and save results
def run_crossmatch(fileA, fileB, out_dir="crossmatch_results", max_trees=50):
    os.makedirs(out_dir, exist_ok=True)

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