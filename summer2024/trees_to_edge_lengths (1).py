import dendropy as dpy
from dendropy.simulate import treesim
import os
import pandas as pd
from subprocess import *


#   Notes: 
#       - should identify non-binary and star tree topologies
#       - assumes gene trees are rooted
#       - the code to generate each output file is only run if that output files does not already exist
#       - variables set by the user for all scenarios (simulating gene trees and not) 
#           are at the top in all caps the user sets are at the top in all caps

# OUTPUT
#       - CSV file with splits as columns, edge lengths as values, one tree per line (same order as Newick gene tree file)
#               ex. gt_edge_lengths_CAT_tauAB-10000_tauABC-20000_tauRoot-30000_N1-10000_N2-10000_N3-10000.csv
#

# directory for output with / at end
OUTDIR = "/Users/kojoa/Downloads/tree-hypothesis-testing-main/tree-hypothesis-testing-main/TreeData/simdata2/"

SIMULATE_GENE_TREES = True      # set to false if using an existing gene tree file
                                # Taxa will be A, B, C, D

CALCULATE_FRECHET_MEAN = False # set to false to avoid calculating the Frechet mean, even if it hasn't already been calculated

# Location of jar file for computing Frechet mean
JAR_FILE = "SturmMean_201102.jar"


################################################
# Parameters for using existing gene tree file
################################################

# Gene tree file
GENE_TREE_FILE = "/Users/megan/Dropbox/RESEARCH/AIM_Square/BHVsummaries/BHV-contours/Gibbons/Bal-noSSY_generations.txt"

# Taxa order will correspond to A, B, C, D, where species tree is either (((A,B),C),D) or ((A,B),(C,D))
TAXA_LIST = ["HPL","HMO","HLE","NLE"]

# Character to put between taxa in the column names
CHAR_BTW_TAXA = "-"


####################################################
# Parameters for simulating new gene tree file using the dendropy package from either a 
#   caterpillar or balanced 4-leaf species tree.
# File name will contain these parameter settings and taxa_list will be set to A, B, C, D
#
#   User inputs (set variables below):
#       - tau1, tau2, tau3  (timings of speciations in species tree; see paper for details)
#       - N1, N2, N3 (population sizes for branches in species tree; see paper for details)
#       - number of gene trees to simulate
#       - whether species tree is caterpillar or balanced   
#
#
#   Output (in user-specified outdir directory):
#       - gene trees, one per line in Newick
#               ex. gts_dendropy_CAT_tauAB-10000_tauABC-20000_tauRoot-30000_N1-10000_N2-10000_N3-10000.tre
#
#   Returns: name of the gene tree file
####################################################

CATERPILLAR_SPECIES_TREE = True  # set to false for balanced species tree

TAU1 = 10000  # tau_AB
TAU2 = 20000  # tau_ABC
TAU3 = 30000   # tau_root

POP_SIZE = 10000  # default population size 
N1 = 10000    # population size above AB
N2 = 10000     # population size above ABC
N3= 15000     # population size above root

NUM_GENE_TREES = 10000


#####################################################
#####################################################
#
#  End of user settings
#
#####################################################
#####################################################



#############################################################
#
#  Script
#
############################################################

def main():

    # Simulate the gene trees, if desired
    if (SIMULATE_GENE_TREES):

        if (CATERPILLAR_SPECIES_TREE):
            tau_N_settings = "_tauAB-" + str(TAU1) + "_tauABC-" + str(TAU2) + "_tauRoot-"+ str(TAU3) + "_pAB-" + str(N1) + "_pABC-" + str(N2) + "_pRoot-" + str(N3)
            tree_type = "CAT"

        else: 
            tau_N_settings = "_tauAB-" + str(TAU1) + "_tauCD-" + str(TAU2) + "_tauRoot-"+ str(TAU3) + "_pAB-" + str(N1) + "_pCD-" + str(N2) + "_pRoot-" + str(N3)
            tree_type = "BAL" 

        # Create the gene tree file name based on the parameters
        gene_tree_file = OUTDIR + "gts_dendropy_" + tree_type + tau_N_settings
        
        simulate_gene_trees(TAU1, TAU2, TAU3, N1, N2, N3, POP_SIZE, NUM_GENE_TREES, CATERPILLAR_SPECIES_TREE, gene_tree_file)

        taxa_list = ["A","B","C","D"]
        char_btw_taxa = ""
    else:
        gene_tree_file = GENE_TREE_FILE
        taxa_list = TAXA_LIST
        char_btw_taxa = CHAR_BTW_TAXA
        
    (head, out_suffix) = os.path.split(gene_tree_file)

    csv_file = OUTDIR + "edge_lengths_" + out_suffix + ".csv"


    ####################################################################################
    # Extract Edge Lengths
    #
    # Create a dataframe with each row representing a gene tree,
    # and each column representing a clade and containing the edge length in the tree
    ####################################################################################

    if not os.path.exists(csv_file):
        newick_to_csv(gene_tree_file, csv_file, character_btw_taxa = char_btw_taxa, taxa_list = taxa_list)

    else:
        print("Skipping generating CSV file of edge lengths -",csv_file,"already exists")


    ######################################################################################
    #
    # Compute Frechet Mean and extract edge lengths to CSV, if mean is wanted and not already calculated
    # 
    # Using stopping criteria (same as in paper with Dan Brown):
    #   cauchySeqLen = 10
    #   epsilon = 0.000001
    #   maxIterations = 5000000
    #####################################################################################

    if (CALCULATE_FRECHET_MEAN):
        frechet_mean_file = OUTDIR + "frechet_mean_" + out_suffix
        mean_csv = OUTDIR + "frechet_mean_" + out_suffix + ".csv"

        # Compute the mean
        if not os.path.exists(frechet_mean_file):
    
            compute_frechet_mean(gene_tree_file, frechet_mean_file)

        else:
            print("Skipping computing Frechet mean; either not required or",frechet_mean_file,"already exists")

        # Extract the edge lengths
        if not os.path.exists(mean_csv):

            extract_frechet_mean(frechet_mean_file, mean_csv, character_btw_taxa = char_btw_taxa, taxa_list = taxa_list)

        else: 
            print("Skipping generating the CSV file containing information about the Frechet mean; either not required or",mean_csv,"already exists")



#################################################################
# Simulate gene trees from the species tree using the coalescent
#
# Parameters: 
#       - tau1, tau2, tau3  (timings of speciations in species tree; see paper for details)
#       - N1, N2, N3 (population sizes for branches in species tree; see paper for details)
#       - number of gene trees to simulate
#       - whether species tree is caterpillar or balanced 
#       - output directory for the gene tree file (default is directory containing script)
#################################################################

def simulate_gene_trees(tau1, tau2, tau3, N1, N2, N3, pop_size, num_gene_trees, caterpillar_species_tree, gene_tree_file):
    if caterpillar_species_tree:
        # Caterpillar species tree is (((A,B),C),D)
        # With branch lengths given by taus:
        # (((A:tau1,B:tau1):(tau2 - tau1),C:tau2):(tau3 - tau2),D:tau3);
        sp_tree_str = "[&R] (((A:{},B:{}):{},C:{}):{},D:{});".format(tau1,tau1,tau2 - tau1,tau2,tau3 - tau2,tau3)
        
        # Define the population sizes for the edges, in the same order the edges' lengths are
        # encountered in the Newick string. 
        pop_sizes = [pop_size,          #A      
                    pop_size,       # B
                    N1,             # AB
                    pop_size,       # C
                    N2,             # ABC
                    pop_size,       # D
                    N3]             # ABCD (root)

    else:
        # Balanced species tree is ((A,B),(C,D))
        # With branch lengths given by taus:
        # ((A:tau1,B:tau1):(tau3 - tau1),(C:tau2,D:tau2):(tau3 - tau2));
        sp_tree_str = "[&R] ((A:{},B:{}):{},(C:{},D:{}):{});".format(tau1,tau1,tau3-tau1,tau2,tau2,tau3-tau2)

        # Define the population sizes for the edges, in the same order the edges' lengths are
        # encountered in the Newick string. 
        pop_sizes = [pop_size,          #A      
                    pop_size,       # B
                    N1,             # AB
                    pop_size,       # C
                    pop_size,       # D
                    N2,             # CD
                    N3]             # ABCD (root)

    if os.path.exists(gene_tree_file):
        print("Skipping simulating gene trees - ",gene_tree_file, " already exists")
        return

    print("Species tree is",sp_tree_str)

    # Create the species tree as a dendropy object
    sp_tree = dpy.Tree.get(data = sp_tree_str, schema = "newick")
    
    # Assign population sizes to the edges
    for (edge, pop_size) in zip(sp_tree.postorder_edge_iter(), pop_sizes):
        edge.pop_size = pop_size

    # Create a taxon mapping to use for the gene trees
    gene_to_species_map = dpy.TaxonNamespaceMapping.create_contained_taxon_mapping(
            containing_taxon_namespace=sp_tree.taxon_namespace,
            num_contained=[1, 1, 1, 1])

    # Simulate gene trees, link to the species tree taxon mapping, and write the output file
    with open(gene_tree_file,"w") as fout:

        for i in range(num_gene_trees):
            gene_tree = treesim.contained_coalescent_tree(containing_tree=sp_tree,
                gene_to_containing_taxon_map=gene_to_species_map, default_pop_size =pop_size)
            
            # Remove the _1 in the gene taxa names
            fout.write(gene_tree.as_string(schema='newick').replace("_1",""))

    print("Finished simulating gene trees in", gene_tree_file)

    return(gene_tree_file)


################################################################
# Read in trees, either from a Newick file (default) or 
# from a string of a single tree in Newick format (pass in False for in_trees_are_file parameter)
# and create a CSV file with the edge lengths.
#
# For rooted trees with 4 leaves only.
#
# Parameters:
# in_trees:  either a file name of tree file or a Newick string of a single tree
# csv_file:  CSV output file name
# character_btw_taxa: character to put between taxa names in CSV column names; - by default
# in_trees_are_file: True if in_trees is a file name; False if in_trees is a Newick string
# taxa_order: list of taxa in the order (ie. alphabetical) they should be used in the CSV column names; 
#                by default the taxa can be ordered arbitrarily
#################################################################

def newick_to_csv(in_trees,csv_file,character_btw_taxa = "-", in_trees_are_file = True, taxa_list = None):

    # If specifying taxa order, make the namespace from the list
    if (taxa_list != None):
        taxon_namespace = dpy.TaxonNamespace(taxa_list, label = "taxa")

        # Read in the tree or trees
        if in_trees_are_file:
            # Read in the tree file
            # Read all gene trees into a Dendropy TreeList (common taxon namespace)
            trees = dpy.TreeList.get(
                path = in_trees,
                taxon_namespace = taxon_namespace,
                schema= "newick",
                rooting= 'default-rooted')
        else:
            # Read in the single tree
            trees = dpy.TreeList.get(
                data = in_trees,
                taxon_namespace = taxon_namespace,
                schema= "newick",
                rooting= 'default-rooted')
    else:
        # Read in the tree or trees and have the namespace generated arbitrarily
        if in_trees_are_file:
            # Read in the tree file
            # Read all gene trees into a Dendropy TreeList (common taxon namespace)
            trees = dpy.TreeList.get(
                path = in_trees,
                schema= "newick",
                rooting= 'default-rooted')
        else:
            # Read in the single tree
            trees = dpy.TreeList.get(
                data = in_trees,
                schema= "newick",
                rooting= 'default-rooted')

        taxa_list = [str(x).replace('`','').replace("'","").replace(" ","_") for x in trees[0].taxon_namespace]

    #print("Taxa list is", taxa_list)

    # Count number of trees
    num_trees = sum(1 for _ in trees)
 
    # Create the dataframe for the edge lengths
    # The four bits in the column names represent the four taxa (as ordered in the dendropy taxa name space), 
    # and 1's indicate which taxa are in that column's clade.
    cols = ["0001", "0010", "0100", "1000","0011", "0101","0110", "1001", "1010", "1100",\
            "0111","1011","1101","1110", "topology"]
    df = pd.DataFrame(0, index= range(0,num_trees), columns = cols)

    # For each tree, get its edge lengths
    num_line = 0
    for tree in trees:
        #print()
        #print("tree is", tree)
        #print("taxa are", tree.taxon_namespace)
        for bitmask,edge in tree.bipartition_edge_map.items():
            #print(bitmask,edge.length)
            # Skip the root
            if (str(bitmask) == '1111'): continue
            df.loc[num_line,str(bitmask)]= edge.length
        num_line +=1

    # Get the topology of the tree in each row
    # The code below, which makes a list of all the topologies and then copies it into the dataframe
    # is faster than placing each topology in the dataframe as we go (ie. df.loc[i,'topology'] = '(((AB)C)D)' , etc.)
    
    tp = []

    # conversion:  1000 -> D -> taxa_list[3]
    #			   0100 -> C -> taxa_list[2]
    #              0010 -> B -> taxa_list[1]
    #              0001 -> A -> taxa_list[0]
    for i in range(0,num_trees):
        if df.loc[i,'0011']!=0:             # AB
            if df.loc[i,'0111']!=0:
                topology = '(((A,B),C),D)'
            elif df.loc[i,'1011']!=0:
                topology = '(((A,B),D),C)'
            elif df.loc[i,'1100']!=0:
                topology = '((A,B),(C,D))'
            else:
                topology = '((A,B),C,D)'
        elif df.loc[i,'0101']!=0:           # AC
            if df.loc[i,'0111']!=0:
                topology = '(((A,C),B),D)'
            elif df.loc[i,'1101']!=0:
                topology = '(((A,C),D),B)'
            elif df.loc[i,'1010'] != 0:
                topology = '((A,C),(B,D))'
            else:
                topology = '((A,C),B,D)'
        elif df.loc[i,'1001']!=0:           # AD
            if df.loc[i,'1011']!=0:
                topology = '(((A,D),B),C)'
            elif df.loc[i,'1101']!=0:
                topology = '(((A,D),C),B)'
            elif df.loc[i,'0110']!=0:
                topology = '((A,D),(B,C))'
            else:
                topology = '((A,D),B,C)'
        elif df.loc[i,'0110']!=0:           # BC
            if df.loc[i,'0111']!=0:
                topology = '(((B,C),A),D)'
            elif df.loc[i,'1110']!=0:
                topology = '(((B,C),D),A)'
            else:				
                topology = '((B,C),A,D)'
        elif df.loc[i,'1010']!=0:           # BD
            if df.loc[i,'1011'] != 0:
                topology = '(((B,D),A),C)'
            elif df.loc[i,'1110']:
                topology = '(((B,D),C),A)'
            else:
                topology = '((B,D),A,C)'
        elif df.loc[i,'1100']!=0:           # CD
            if df.loc[i,'1101']!=0:
                topology = '(((C,D),A),B)'
            elif df.loc[i,'1110'] != 0:
                topology = '(((C,D),B),A)'
            else:
                topology = '((C,D),A,B)'
        elif df.loc[i,'0111'] !=0:          # ABC
            topology = '((A,B,C),D)'
        elif df.loc[i,'1011'] !=0:          # ABD
            topology = '((A,B,D),C)'
        elif df.loc[i,'1101'] != 0:         # ACD
            topology = '((A,C,D),B)'
        elif df.loc[i,'1110'] != 0:         # BCD
            topology = '((B,C,D),A)'
        else:
            topology = '(A,B,C,D)'
        i += 1
        taxa_dict = {'A': taxa_list[0], 
                    'B':taxa_list[1], 
                    'C': taxa_list[2],
                    'D': taxa_list[3]}
        topology = topology.translate(str.maketrans(taxa_dict))
        tp.append(topology)
    df['topology'] = tp

    # Convert all bitmap column names into ones with letters, according to the trees' common namespace
    new_col_names = []
    for col in cols:
        if col == "topology": continue
        label = []
        for idx,bit in enumerate(reversed(col)):
            if bit=='1':
                label.append(taxa_list[idx])
        new_col_names.append(character_btw_taxa.join(label))
    new_col_names.append("topology")

    df.columns = new_col_names

    df.to_csv(csv_file,index = False)

######################################################################################
#
# Compute Frechet Mean
# 
# Using stopping criteria (same as in paper with Dan Brown):
#   cauchySeqLen = 10
#   epsilon = 0.000001
#   maxIterations = 5000000
#####################################################################################

def compute_frechet_mean(tree_file, frechet_mean_file):
    command = "java -jar " + JAR_FILE + " -a rand_perm -o " + frechet_mean_file + " -c 10 -e 0.000001 -n 5000000 " + tree_file
    print("Calling: " + command)
    call(command.split())

########################################################################
#
#  Extract Frechet Mean
#
# Need to extract the Newick representation for the Frechet mean from the output file
# generated by computing the Frechet mean
# and store the edge length information in a CSV file.
#########################################################################

def extract_frechet_mean(frechet_mean_file, mean_csv, character_btw_taxa = "-", taxa_list = None):
    with open(frechet_mean_file,'r') as fin:
        all_lines = fin.readlines()
    mean_tree = all_lines[4]   # get the fifth line in the Frechet mean file

    newick_to_csv(mean_tree,mean_csv,in_trees_are_file = False, character_btw_taxa = character_btw_taxa, taxa_list = taxa_list)


###########################################################################
#
#  Scale tree edge lengths
#
#   Multply all tree edge lengths in all trees in the file by x.
#
###########################################################################

def scale_trees(in_file, out_file, x):
    trees = dendropy.Tree.get(
                path= in_file,
                schema='newick',
                rooting= 'default-rooted')

    for tree in trees:
        for edge in tree.postorder_edge_iter():
            if edge.length is None:
                edge.length = 0
            else:
                edge.length = float(edge.length)*x

    trees.write(path = out_file,
                schema = "newick")


if __name__ == "__main__":
    main()
    