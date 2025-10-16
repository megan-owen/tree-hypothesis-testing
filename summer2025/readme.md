- SturmMean_201102.jar
Java tool for computing Fréchet (Sturm) means of a set of trees.

- gtp.jar
Java tool for BHV tree-space geodesics (used via JPype or CLI to compute BHV distances).

- trees_to_edge_lengths.py
Simulate gene trees (DendroPy), compute edge lengths via SturmMean (Fréchet mean workflow), and write outputs.

- mean_test_gene_tree.py
Difference-of-means permutation test on edge-length features derived from the Fréchet means.

- crossmatch_functions.py
Implementation of the cross-match test on a set of items given

- crossmatch_test_gene_trees.py
Cross-match test using weighted Robinson–Foulds (wRF) distances from DendroPy.

 - crossmatch_bhv.py
Cross-match test using BHV distances computed via gtp.jar through JPype.

- BHV FILE
Practice Calls using Jpype calling Java classes we already have in a python script.

Tested Crossmatches for Locus trees using BHV and Weighted Robinsons-Foulds but decided to stick with Robinsons-Fould since there wasn't much difference in how much time it took to compute them both.




- The mean_test_results/folder contains output from the permutation (mean) test performed across 304 loci.

Graphs: Histograms are available for loci 1–75 and 152–304.

CSV Results: The summary CSV (mean_test_summary.csv) includes data for all loci (1–304).
These results summarize the p-values and test statistics used to compare the two gene tree sets.



HOW TO GUIDE FOR GENE TREES:

To generate two sets of gene trees using trees_to_edge_lengths.py (which uses argparse), you can run two separate commands from the command line(terminal) like this:

python trees_to_edge_lengths.py --tau1 10000 --tau2 10500 --tau3 11000 --N1 10000 --N2 10000 --N3 10000 --out output/

python trees_to_edge_lengths.py --tau1 10000 --tau2 10100 --tau3 11100 --N1 10000 --N2 10000 --N3 10000 --out output/

Output Files Will Be Named Like:

output/gts_dendropy_CAT_tauAB-10000_tauABC-10500_tauRoot-11000_pAB-10000_pABC-10000_pRoot-10000

output/gts_dendropy_CAT_tauAB-10000_tauABC-10100_tauRoot-11100_pAB-10000_pABC-10000_pRoot-10000

Now to run Mean (Permutation testing) put this in the terminal.

python mean_test_gene_tree.py

File names are hardcoded in the mean_test_gene_tree so change them if you dont use the same parametes to generate the gene trees.

mean_test_gene_tree.py output is saved in "mean_test_output.txt"

TO RUN CROSSMATCH TEST WITH WEIGHTED ROBINSON FOULD:

python crossmatch_test_gene_trees.py

Files are also hardcoded in this script. 

TO RUN CROSSMATCH TEST USING BHV RUN:

python crossmatch_bhv.py  

Files are also hardcoded in this script. 










