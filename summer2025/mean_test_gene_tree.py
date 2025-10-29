# Permutation test on Fréchet means:
#   - Loads two samples of gene trees (one Newick per line) from FILE_A / FILE_B.
#   - Uses external helpers to compute Fréchet means per group and a weighted RF
#     distance between those means as the test statistic.
#   - Shuffles trees across groups to build a null distribution; p-value is
#     Pr[ permuted distance >= observed distance ].
#
# Notes:
#   * This script uses a file-based Fréchet-mean helper:
#       - compute_frechet_mean(input_tree_file, output_mean_file)
#       - extract_frechet_mean(output_mean_file) -> Newick string
#   * The distance between means is computed via weighted RF (eric_functions.weighted_distance).
import random
import shutil
import matplotlib.pyplot as plt
import os

from eric_functions import weighted_distance
from trees_to_edge_lengths import compute_frechet_mean, extract_frechet_mean

# =============================================================================
#                         USER-CONFIGURABLE VARIABLES
# =============================================================================

# Input tree files (one Newick tree per line)
FILE_A = "output/gts_dendropy_CAT_tauAB-10000.0_tauABC-10500.0_tauRoot-11000.0_pAB-10000_pABC-10000_pRoot-10000"
FILE_B = "output/gts_dendropy_CAT_tauAB-10000.0_tauABC-10100.0_tauRoot-11100.0_pAB-10000_pABC-10000_pRoot-10000"

# Precomputed Fréchet means for each sample
MEAN_A = "output/frechet_mean_gts_dendropy_CAT_tauAB-10000.0_tauABC-10500.0_tauRoot-11000.0_pAB-10000_pABC-10000_pRoot-10000"
MEAN_B = "output/frechet_mean_gts_dendropy_CAT_tauAB-10000.0_tauABC-10100.0_tauRoot-11100.0_pAB-10000_pABC-10000_pRoot-10000"

# Number of permutations and output folder
NUM_PERMUTATIONS = 100
SAVE_PATH = "new_mean_test_results"


# =============================================================================
#                                MAIN ROUTINE
# =============================================================================


# Run a permutation test using the distance between Fréchet means.
# Parameters:
#   fileA (str): path to sample A file (one Newick per line).
#   fileB (str): path to sample B file (one Newick per line).
#   real_mean_fileA (str): path to precomputed Fréchet mean for A.
#   real_mean_fileB (str): path to precomputed Fréchet mean for B.
#   save_path (str): output directory (stores temp files and final report).
#   jar_file (str): unused here; kept for API parity with other pipelines.
#   num_permutations (int): number of random shuffles for the null distribution.
#   plot (bool): if True, show a histogram of permuted distances with observed marked.
# Returns:
#   dict: {"observed_diff": float, "p_value": float}

def run_mean_test(fileA, fileB, real_mean_fileA, real_mean_fileB,save_path="new_mean_test_results",
    jar_file="SturmMean_201102.jar",num_permutations=100, plot=False):

    os.makedirs(save_path, exist_ok=True)
    temp_dir = os.path.join(save_path, "tmp_perm")
    os.makedirs(temp_dir, exist_ok=True)

    #Temporary files used per permutation
    tree_file1 = os.path.join(temp_dir, "group1.tre")
    tree_file2 = os.path.join(temp_dir, "group2.tre")
    mean_file1 = os.path.join(temp_dir, "mean1.tre")
    mean_file2 = os.path.join(temp_dir, "mean2.tre")

    # Load tree sammples from input files 
    with open(fileA, 'r') as f:
        trees1 = [line.strip() for line in f if line.strip()]
    with open(fileB, 'r') as f:
        trees2 = [line.strip() for line in f if line.strip()]

    all_trees = trees1 + trees2
    labels = [0] * len(trees1) + [1] * len(trees2)

    # Compute observed test statistic (distance between the true frechet means) from real Frechet means 
    mean1 = extract_frechet_mean(real_mean_fileA)
    mean2 = extract_frechet_mean(real_mean_fileB)
    observed_distance = weighted_distance(mean1, mean2)

    distances = []

    for i in range(num_permutations):
        print(f"\n=== Running permutation {i+1}/{num_permutations} ===")

        # Shuffle trees and split into two groups 
        shuffled = list(zip(all_trees, labels))
        random.shuffle(shuffled)
        group1 = [tree for tree, _ in shuffled[:len(trees1)]]
        group2 = [tree for tree, _ in shuffled[len(trees1):]]

        # Write the shuffled groups to temporary files
        with open(tree_file1, 'w') as f:
            f.write("\n".join(group1) + "\n")
        with open(tree_file2, 'w') as f:
            f.write("\n".join(group2) + "\n")

        # Compute Frechet means and distance for the permuted groups
        compute_frechet_mean(tree_file1, mean_file1)
        compute_frechet_mean(tree_file2, mean_file2)

        perm_mean1 = extract_frechet_mean(mean_file1)
        perm_mean2 = extract_frechet_mean(mean_file2)
        dist = weighted_distance(perm_mean1, perm_mean2)
        distances.append(dist)

    #compute the p-value 
    p_value = sum(d >= observed_distance for d in distances) / num_permutations

    # Save results 
    result_file = os.path.join(save_path, "new_mean_test_output.txt")
    with open(result_file, "w") as f:
        f.write("Permutation Test on Frechet Means:\n")
        f.write("----------------------------------\n")
        f.write(f"Observed Frechet mean distance: {observed_distance:.6f}\n")
        f.write(f"P-value from permutation test: {p_value:.6f}\n")
        if p_value < 0.05:
            f.write("Null hypothesis rejected at 0.05 level.\n")
        else:
            f.write("Failed to reject the null hypothesis at 0.05 level.\n")

    if plot:
        plt.hist(distances, bins=30, color="skyblue", edgecolor="black")
        plt.axvline(observed_distance, color="red", linestyle="--", label=f"Observed = {observed_distance:.2f}")
        plt.title("Permutation Test: Distance Between Frechet Means")
        plt.xlabel("Distance")
        plt.ylabel("Frequency")
        plt.legend()
        plt.tight_layout()
        plt.show()

    # Clean up and delete temporary files
    shutil.rmtree(temp_dir)

    return {
        "observed_diff": observed_distance,
        "p_value": p_value
    }

# === Run Test ===
if __name__ == "__main__":
    run_mean_test(FILE_A, FILE_B, MEAN_A, MEAN_B, plot=True)

