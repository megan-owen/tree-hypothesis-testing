"""
mean_test_gene_tree.py

This script defines the `run_mean_test` function, which performs a permutation test
to assess whether two sets of gene trees differ in their average edge lengths.

Function: run_mean_test(fileA, fileB, save_path="mean_test_results", plot=False)
- Inputs:
    - fileA, fileB: Paths to CSV files containing edge lengths from two tree samples.
    - save_path: Directory where results will be saved.
    - plot: If True, saves a histogram of the null distribution.
- Method:
    - Computes observed difference in mean edge lengths.
    - Performs 10,000 permutations to build a null distribution.
    - Calculates a p-value by comparing the observed difference to the null.
- Outputs:
    - Text summary saved in mean_test_results/mean_test_output.txt.
    - (Optional) Histogram saved as histogram_mean_diff.png.

Returns:
- A dictionary with:
    - "observed_diff": observed difference in means
    - "p_value": p-value from the permutation test
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

def run_mean_test(fileA, fileB, save_path="mean_test_results", plot=False):
    os.makedirs(save_path, exist_ok=True)

    dfA = pd.read_csv(fileA).drop(columns=["topology"]).astype(float)
    dfB = pd.read_csv(fileB).drop(columns=["topology"]).astype(float)

    # Compute the observed mean difference
    meanA = dfA.mean().mean()
    meanB = dfB.mean().mean()
    observed_diff = abs(meanA - meanB)

    # Combine both datasets and prepare for permutation test
    combined = pd.concat([dfA, dfB], ignore_index=True)
    n_A = len(dfA)
    num_simulations = 10000
    simulated_diffs = []

    for i in range(num_simulations):
        shuffled = combined.sample(frac=1, random_state=i).reset_index(drop=True)
        groupA_sim = shuffled.iloc[:n_A]
        groupB_sim = shuffled.iloc[n_A:]
        sim_diff = abs(groupA_sim.mean().mean() - groupB_sim.mean().mean())
        simulated_diffs.append(sim_diff)

    p_value = np.mean(np.array(simulated_diffs) >= observed_diff)

    # Plot the histogram if requested
    if plot:
        plt.hist(simulated_diffs, bins=50, alpha=0.7, label="Null distribution")
        plt.axvline(observed_diff, color="red", linestyle="--", label="Observed diff")
        plt.xlabel("Difference in Means")
        plt.ylabel("Frequency")
        plt.title("Permutation Test on Edge Lengths of Gene Trees")
        plt.legend()
        plt.tight_layout()
        plt.savefig(f"{save_path}/histogram_mean_diff.png")
        plt.close()

    # Save numerical result
    with open(f"{save_path}/mean_test_output.txt", "w") as f:
        f.write("Difference of Means Test:\n")
        f.write("-------------------------\n")
        f.write(f"Observed mean difference: {observed_diff:.4f}\n")
        f.write(f"P-value from permutation test: {p_value:.4f}\n")
        if p_value < 0.05:
            f.write("Null hypothesis rejected at 0.05 level.\n")
        else:
            f.write("Failed to reject the null hypothesis at 0.05 level.\n")

    return {
        "observed_diff": observed_diff,
        "p_value": p_value
    }