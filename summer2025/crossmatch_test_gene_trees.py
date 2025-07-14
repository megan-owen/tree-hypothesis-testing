import dendropy
import numpy as np
import os
from crossmatch_functions import crossmatchtest
from dendropy.calculate import treecompare  # import for RF distances

# Load trees from a .tre file
def load_trees(filename):
    return dendropy.TreeList.get(path=filename, schema="newick", preserve_underscores=True)

#Temporarily testing with only 50 trees from each sample.
sampleA_trees = load_trees("sampleA.tre")[:50]
sampleB_trees = load_trees("sampleB.tre")[:50]


# Combine both samples into one list
all_trees = sampleA_trees + sampleB_trees
num_A = len(sampleA_trees)
num_B = len(sampleB_trees)


print(f"Loaded {num_A} trees from sampleA and {num_B} from sampleB.")
print(f"Total trees (all_trees): {len(sampleA_trees + sampleB_trees)}")

# Create label vector: 0 for sample A, 1 for sample B
labels = [0] * num_A + [1] * num_B

# Initialize pairwise distance matrix
n = num_A + num_B
distance_matrix = np.zeros((n, n))

# Compute pairwise unweighted Robinson-Foulds distances
comparison_count = 0
total_comparisons = (len(all_trees) * (len(all_trees) - 1)) // 2
print(f"Starting pairwise RF distance calculations ({total_comparisons} comparisons)...")

for i in range(len(all_trees)):
    for j in range(i + 1, len(all_trees)):
        dist = treecompare.symmetric_difference(all_trees[i], all_trees[j])
        distance_matrix[i][j] = dist
        distance_matrix[j][i] = dist

        comparison_count += 1
        if comparison_count % 1000 == 0:
            print(f"Completed {comparison_count}/{total_comparisons} comparisons...")

# Run the cross match statistical test
raw_result = crossmatchtest(labels, distance_matrix)
 
result = {
    "a1": raw_result[0],
    "Ea1": raw_result[1],
    "Va1": raw_result[2],
    "dev": raw_result[3],
    "pval (exact)": raw_result[4],
    "pval (normal approx)": raw_result[5],
}

# Print and save output
print("\nCross-match test results:")
for key, value in result.items():
    print(f"{key}: {value}")

os.makedirs("crossmatch_results", exist_ok=True)
with open("crossmatch_results/crossmatch_output.txt", "w") as f:
    for key, value in result.items():
        f.write(f"{key}: {value}\n")
