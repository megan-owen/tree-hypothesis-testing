import os
import time
import numpy as np
import jpype
import jpype.imports
from crossmatch_functions import crossmatchtest

start_time = time.time()

if not jpype.isJVMStarted():
    # Make sure gtp.jar is in the working dir or give an absolute path
    jpype.startJVM(classpath=["gtp.jar"])

PolyMain = jpype.JClass("polyAlg.PolyMain")
PhyloTree = jpype.JClass("distanceAlg1.PhyloTree")

# === Inputs ===
FILE_A = "output/gts_dendropy_CAT_tauAB-10000.0_tauABC-10500.0_tauRoot-11000.0_pAB-10000_pABC-10000_pRoot-10000"
FILE_B = "output/gts_dendropy_CAT_tauAB-10000.0_tauABC-10100.0_tauRoot-11100.0_pAB-10000_pABC-10000_pRoot-10000"
MAX_TREES = 300
OUTPUT_DIR = "crossmatch_results"

# Simple file loader: assumes proper Newick formatting 
def read_newick_file(path, max_trees=None):
    with open(path, "r") as f:
        trees = [line.strip() for line in f if line.strip()]
    return trees[:max_trees] if max_trees else trees

# Compute BHV pairwise distance matrix using Java classes
def compute_bhv_distance_matrix(newick_strings, rooted=True):
    n = len(newick_strings)
    M = np.zeros((n, n), dtype=float)

    print("Computing pairwise BHV distances...")
    for i in range(n):
        for j in range(i + 1, n):
            try:
                t1 = PhyloTree(newick_strings[i], rooted)
                t2 = PhyloTree(newick_strings[j], rooted)
                dist = PolyMain.getGeodesic(t1, t2, None).getDist()
                M[i, j] = M[j, i] = dist
            except Exception as e:
                print(f"Error computing distance between {i} and {j}: {e}")
                M[i, j] = M[j, i] = np.nan
        if (i + 1) % 25 == 0 or i == n - 1:
            print(f"  row {i+1}/{n} done")
    return M

# Run crossmatch test
def run_crossmatch(fileA, fileB, out_dir="crossmatch_results", max_trees=50, rooted=True):
    os.makedirs(out_dir, exist_ok=True)

    treesA = read_newick_file(fileA, max_trees)
    treesB = read_newick_file(fileB, max_trees)

    print(f"Loaded {len(treesA)} trees from {fileA}")
    print(f"Loaded {len(treesB)} trees from {fileB}")

    all_trees = treesA + treesB
    labels = [0] * len(treesA) + [1] * len(treesB)

    print("Computing BHV distance matrix...")
    M = compute_bhv_distance_matrix(all_trees, rooted=rooted)

    print("Running cross-match test...")
    a1, Ea1, Va1, dev, pval_exact, pval_normal = crossmatchtest(labels, M)

    results = {
        "a1": a1,
        "Ea1": Ea1,
        "Va1": Va1,
        "dev": dev,
        "pval (exact)": pval_exact,
        "pval (normal approx)": pval_normal,
    }

    base1 = os.path.basename(fileA)
    base2 = os.path.basename(fileB)
    out_file = os.path.join(out_dir, f"crossmatch_BHV_{base1}_vs_{base2}.txt")

    with open(out_file, "w") as f:
        for k, v in results.items():
            print(f"{k}: {v}")
            f.write(f"{k}: {v}\n")

    print(f"\nResults saved to: {out_file}")
    return results

# --- Main ---
if __name__ == "__main__":
    run_crossmatch(FILE_A, FILE_B, out_dir=OUTPUT_DIR, max_trees=MAX_TREES, rooted=True)
    print(f"\nTotal time: {time.time() - start_time:.2f} seconds")
    jpype.shutdownJVM()
