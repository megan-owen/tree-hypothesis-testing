import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load the full edge length data (gene tree samples)
fileA = "edge_lengths_gts_dendropy_CAT_tauAB-10000.0_tauABC-10100.0_tauRoot-11100.0_pAB-10000_pABC-10000_pRoot-10000.csv"
fileB = "edge_lengths_gts_dendropy_CAT_tauAB-10000.0_tauABC-10100.0_tauRoot-11100.0_pAB-20000_pABC-20000_pRoot-20000.csv"

# Read the CSVs
dfA = pd.read_csv(fileA).drop(columns=["topology"])
dfB = pd.read_csv(fileB).drop(columns=["topology"])

# Convert to floats to avoid warnings
dfA = dfA.astype(float)
dfB = dfB.astype(float)

# Calculate the observed difference in means
meanA = dfA.mean().mean()
meanB = dfB.mean().mean()
observed_diff = abs(meanA - meanB)

print(f"\nObserved mean difference in edge lengths: {observed_diff:.4f}")

#Run a permutation test 
combined = pd.concat([dfA, dfB], ignore_index=True)

n_A = len(dfA)
n_B = len(dfB)

num_simulations = 10000
simulated_diffs = []

for i in range(num_simulations):
    shuffled = combined.sample(frac=1, random_state=i).reset_index(drop=True)
    groupA_sim = shuffled.iloc[:n_A]
    groupB_sim = shuffled.iloc[n_A:]
    
    sim_diff = abs(groupA_sim.mean().mean() - groupB_sim.mean().mean())
    simulated_diffs.append(sim_diff)

#  Calculate p-value
#A low p-value ( < 0.05) means you reject the null. The two sets of gene trees likely come from different distributions.
p_value = np.mean(np.array(simulated_diffs) >= observed_diff)
print(f"P-value from permutation test: {p_value:.4f}\n")



# Plot the results
plt.hist(simulated_diffs, bins=50, alpha=0.7, label="Null distribution")
plt.axvline(observed_diff, color="red", linestyle="--", label="Observed diff")
plt.xlabel("Difference in Means")
plt.ylabel("Frequency")
plt.title("Permutation Test on Edge Lengths of Gene Trees")
plt.legend()
plt.tight_layout()
plt.show()


# Save histogram plot
plt.savefig("mean_test_results/histogram_mean_diff.png")  # Saves the histogram image

# Save the numeric result (test statistic and p-value) to a text file
with open("mean_test_results/mean_test_output.txt", "w") as f:
    f.write(f"Observed mean difference in edge lengths: {observed_diff:.4f}\n")
    f.write(f"P-value from permutation test: {p_value:.4f}\n")

