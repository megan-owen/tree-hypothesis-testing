# Crossmatch two-sample test over gene trees using a BHV (Billera–Holmes–Vogtmann)
# distance matrix computed via gtp.jar (Java) accessed through JPype.
#
# Workflow:
#   1) Load two Newick sets (FILE_A, FILE_B; one tree per line).
#   2) Start a JVM with gtp.jar on the classpath and bind PolyMain / PhyloTree.
#   3) Build a symmetric pairwise BHV geodesic distance matrix on A∪B.
#   4) Run the Crossmatch test and write a human-readable report.

import os
import time
import numpy as np 
import jpype
import jpype.imports
from crossmatch_functions import crossmatchtest

# Track total runtime
start_time = time.time()

# =============================================================================
#                           JVM INITIALIZATION
# =============================================================================
# Start the JVM once with gtp.jar on the classpath (required for BHV distances).
if not jpype.isJVMStarted():
    # Make sure gtp.jar is in the working dir or give an absolute path
    jpype.startJVM(classpath=["gtp.jar"])

# Bind Java classes for geodesic computations
PolyMain = jpype.JClass("polyAlg.PolyMain")
PhyloTree = jpype.JClass("distanceAlg1.PhyloTree")

# =============================================================================
#                         USER-CONFIGURABLE VARIABLES
# =============================================================================

# Input tree files (each contains one Newick tree per line)
FILE_A = "output/gts_dendropy_CAT_tauAB-10000.0_tauABC-10500.0_tauRoot-11000.0_pAB-10000_pABC-10000_pRoot-10000"
FILE_B = "output/gts_dendropy_CAT_tauAB-10000.0_tauABC-10100.0_tauRoot-11100.0_pAB-10000_pABC-10000_pRoot-10000"

# Cap the number of trees loaded per file (None = all)
MAX_TREES = 150

# Directory where results are saved
OUTPUT_DIR = "crossmatch_results"

# =============================================================================
#                             HELPER FUNCTIONS
# =============================================================================

# Read Newick trees from file.
# Parameters:
#   path (str): path to .tre or .txt file with one tree per line.
#   max_trees (int|None): maximum trees to load; None = all.
# Returns:
#   list[str]: Newick-formatted tree strings.
def read_newick_file(path, max_trees=None):
    with open(path, "r") as f:
        trees = [line.strip() for line in f if line.strip()]
    return trees[:max_trees] if max_trees else trees


# Compute symmetric pairwise BHV distance matrix using PolyMain.getGeodesic.
# Parameters:
#   newick_strings (list[str]): list of all trees (A + B).
#   rooted (bool): True if trees are rooted (BHV depends on topology type).
# Returns:
#   np.ndarray: (n x n) symmetric matrix of BHV distances.
def compute_bhv_distance_matrix(newick_strings, rooted=True):
    n = len(newick_strings)
    M = np.zeros((n, n), dtype=float)

    print("Computing pairwise BHV distances...")
    for i in range(n):
        for j in range(i + 1, n):
            try:
                # Create Java tree objects and compute BHV geodesic distance
                t1 = PhyloTree(newick_strings[i], rooted)
                t2 = PhyloTree(newick_strings[j], rooted)
                dist = PolyMain.getGeodesic(t1, t2, None).getDist()
                M[i, j] = M[j, i] = dist
            except Exception as e:
                print(f"Error computing distance between {i} and {j}: {e}")
                M[i, j] = M[j, i] = np.nan

        # Progress feedback every 25 trees (or last)
        if (i + 1) % 25 == 0 or i == n - 1:
            print(f"  row {i+1}/{n} done")

    return M


# =============================================================================
#                             CROSSMATCH TEST RUNNER
# =============================================================================

# Run the Crossmatch test on two sets of trees using BHV distances.
# Parameters:
#   fileA (str): path to file with trees for sample A.
#   fileB (str): path to file with trees for sample B.
#   out_dir (str): output directory for result files.
#   max_trees (int): number of trees to load from each file.
#   rooted (bool): whether to treat trees as rooted in BHV computations.
# Returns:
#   dict: crossmatch statistics and p-values.
def run_crossmatch(fileA, fileB, out_dir="crossmatch_results", max_trees=50, rooted=True):
    os.makedirs(out_dir, exist_ok=True)

    # Load samples
    treesA = read_newick_file(fileA, max_trees)
    treesB = read_newick_file(fileB, max_trees)

    print(f"Loaded {len(treesA)} trees from {fileA}")
    print(f"Loaded {len(treesB)} trees from {fileB}")

    # Merge all trees and create label vector (0 for A, 1 for B)
    all_trees = treesA + treesB
    labels = [0] * len(treesA) + [1] * len(treesB)

    # Compute distance matrix
    print("Computing BHV distance matrix...")
    M = compute_bhv_distance_matrix(all_trees, rooted=rooted)

    # Run crossmatch test
    print("Running cross-match test...")
    a1, Ea1, Va1, dev, pval_exact, pval_normal = crossmatchtest(labels, M)

    # Summarize results
    results = {
        "a1": a1,
        "Ea1": Ea1,
        "Va1": Va1,
        "dev": dev,
        "pval (exact)": pval_exact,
        "pval (normal approx)": pval_normal,
    }

    # Save outputs
    base1 = os.path.basename(fileA)
    base2 = os.path.basename(fileB)
    out_file = os.path.join(out_dir, f"crossmatch_BHV_{base1}_vs_{base2}.txt")

    with open(out_file, "w") as f:
        for k, v in results.items():
            print(f"{k}: {v}")
            f.write(f"{k}: {v}\n")

    print(f"\nResults saved to: {out_file}")
    return results


# =============================================================================
#                                 ENTRY POINT
# =============================================================================
if __name__ == "__main__":
    run_crossmatch(FILE_A, FILE_B, out_dir=OUTPUT_DIR, max_trees=MAX_TREES, rooted=True)
    print(f"\nTotal time: {time.time() - start_time:.2f} seconds")

    # Clean JVM shutdown
    jpype.shutdownJVM()
