- SturmMean : 
java jar file that contains code for computing (Frechet) means of a set of trees

- trees_to_edge_lengths.py:

Simulates gene trees and computes edge lengths using the SturmMean Java tool. Produces .tre, .csv, and mean tree summaries.

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






HOW TO GUIDE:

To generate two sets of gene trees using trees_to_edge_lengths.py (which uses argparse), you can run two separate commands from the command line(terminal) like this:

python trees_to_edge_lengths.py --tau1 10000 --tau2 10500 --tau3 11000 --N1 10000 --N2 10000 --N3 10000 --out output/

python trees_to_edge_lengths.py --tau1 10000 --tau2 10100 --tau3 11100 --N1 10000 --N2 10000 --N3 10000 --out output/

Output Files Will Be Named Like:

output/gts_dendropy_CAT_tauAB-10000_tauABC-10500_tauRoot-11000_pAB-10000_pABC-10000_pRoot-10000

output/gts_dendropy_CAT_tauAB-10000_tauABC-10100_tauRoot-11100_pAB-10000_pABC-10000_pRoot-10000



