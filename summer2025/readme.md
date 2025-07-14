- SturmMean : 
java jar file that contains code for computing (Frechet) means of a set of trees

- trees_to_edge_lengths.py:

Simulates gene trees and computes edge lengths using the SturmMean Java tool. Produces .tre, .csv, and mean tree summaries.

- sampleA.tre, sampleB.tre:
Two sets of 10,000 simulated gene trees generated from different population sizes (N=10000 and N=20000 respectively).

- edge_lengths_*.csv:
CSV files containing the edge length vectors of all gene CSV files storing edge lengths for each gene tree sample. Used for the difference-in-means test

- frechet_mean_*.csv:
CSV files summarizing the Fréchet mean tree computed from each gene tree sample.

- mean_test_gene_tree.py:
Implements a difference of means permutation test on edge lengths. Outputs histogram and p-value.
Saves results in mean_test_results/.
updated to save results and reject/fail to reject the null.

- crossmatch_functions.py:
Contains code for the cross-match test based on pairwise RF distances between trees.

- crossmatch_test_gene_trees.p:
Loads samples, computes RF distances, runs the cross-match test, and prints results.

- compare_tests_pipeline.py
(WORK IN PROGRESS) Script to compare the rejection power of both tests by running them 25 times on resampled data. Not yet finalized.