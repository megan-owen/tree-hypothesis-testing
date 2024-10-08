import dendropy
from dendropy.calculate import treecompare
from dendropy.simulate import treesim
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import os
import re
import numpy as np
import random

#folder directory
folder_dir = "/"

def unweighted_distance(tree1,tree2):

    # establish common taxon namespace
    tns = dendropy.TaxonNamespace()

    # ensure all trees loaded use common namespace
    tree1d = dendropy.Tree.get(
            data=tree1,
            schema='newick',
            taxon_namespace=tns)
    tree2d = dendropy.Tree.get(
            data=tree2,
            schema='newick',
            taxon_namespace=tns)

    ## Unweighted Robinson-Foulds distance
    return treecompare.symmetric_difference(tree1d, tree2d)

def weighted_distance(tree1,tree2):

    # establish common taxon namespace
    tns = dendropy.TaxonNamespace()

    # ensure all trees loaded use common namespace
    tree1d = dendropy.Tree.get(
            data=tree1,
            schema='newick',
            taxon_namespace=tns)
    tree2d = dendropy.Tree.get(
            data=tree2,
            schema='newick',
            taxon_namespace=tns)
    return treecompare.weighted_robinson_foulds_distance(tree1d, tree2d)

def euclidean_distance(tree1,tree2):

    # establish common taxon namespace
    tns = dendropy.TaxonNamespace()

    # ensure all trees loaded use common namespace
    tree1d = dendropy.Tree.get(
            data=tree1,
            schema='newick',
            taxon_namespace=tns)
    tree2d = dendropy.Tree.get(
            data=tree2,
            schema='newick',
            taxon_namespace=tns)
    return treecompare.euclidean_distance(tree1d, tree2d)


def generate_sim_tree_files(r_tree2,tree_path):
    sp_tree_str = """\
    [&R] {}""".format(r_tree2)
    
    #loop to generate 304 simulated gene trees
    for i in range(0,304,1):
        
        sp_tree = dendropy.Tree.get(data=sp_tree_str, schema="newick")
        
        #taxon mapping for gene trees to the species tree
        gene_to_species_map = dendropy.TaxonNamespaceMapping.create_contained_taxon_mapping(
                containing_taxon_namespace=sp_tree.taxon_namespace,
                num_contained=1)
        
        gene_tree = treesim.contained_coalescent_tree(containing_tree=sp_tree,
            gene_to_containing_taxon_map=gene_to_species_map,default_pop_size = .06)

        gene_tree_string = gene_tree.as_string(schema='newick')

        #print(gene_tree_string[5:len(gene_tree_string)])

        new_sim_folder = folder_dir + tree_path
        new_sim_file = new_sim_folder + "/SimMltree_{}.txt".format(i+1)
#         if not os.path.exists(new_sim_folder):
#             os.mkdir(folder_dir + "SimMlTrees304_{}".format(copy_num))
        if os.path.exists(new_sim_file):
            os.remove(new_sim_file)
        f= open(new_sim_file,"w+")
        f.write(gene_tree_string[5:len(gene_tree_string)])
        f.close()



def generate_comparison_files(copy_num,sim_file_path,comparison_type):
#     data_dir = "SimMlTrees304_{}/".format(copy_num)

    tree_array = []
    tree_s = ""

    # Load tree strings from files into a list and combine into a single string
    for i in range(1,305,1):
        #change the file pattern here
        file_pat = "/SimMltree_{}.txt".format(i)
        dir_tree = sim_file_path + file_pat 
        tree_string = open(dir_tree).read()
        tree_array.append(tree_string)
        tree_s += tree_string 

    #saving combined string into a new file
    tree_file = open(sim_file_path+ "/SimMltree_merged.txt", "w+")
    tree_file.write(tree_s)
    tree_file.close()    
    
    #CREATE A DATAFRAME TO STORE ALL COMPARISONS OF TREES AND DIFFERENT METHODS OF DISTANCE
    
    col_names = ['Tree_comparison',"Num_First_tree","Num_Second_tree", 'weighted_distance']
    len1 = len(tree_array) - 1
    len2 = int((len1*(len1+1))/2)
    print(len2)
    tree_data = pd.DataFrame(0, index= range(0,len2), columns = col_names)
    #print(tree_data)
    #print("START!!!")
    
    # Loop to compute pairwise distances between the trees
    count = 0
    k = 0
    for i in range(k,304,1):
        for j in range(k+1,304,1):
            #print(str(i+1) + "vs" + str(j+1))
            tree_data.loc[count,'Tree_comparison'] = str(i+1) + "vs" + str(j+1)
            tree_data.loc[count,'Num_First_tree'] = str(i+1)
            tree_data.loc[count,'Num_Second_tree'] = str(j+1)
            
            #for calculating the specified type of tree distance
            if comparison_type == "unweighted_distance":
                tree_data.loc[count,'unweighted_distance'] = unweighted_distance(tree_array[i],tree_array[j])
            elif comparison_type == "weighted_distance":
                tree_data.loc[count,'weighted_distance'] = weighted_distance(tree_array[i],tree_array[j])
            else:
                tree_data.loc[count,'euclidean_distance'] = euclidean_distance(tree_array[i],tree_array[j])
                       
            count += 1    
        k += 1;
    tree_data.to_csv(folder_dir+ 'Sim_Trees_distances_comparisons_{}.csv'.format(copy_num),index=False)

    return tree_data
    #print("END!!")
    
    
    #merged_file : string that contains all the trees.
    
def generate_comparison_files_simphy(merged_file,copy_num,comparison_type):
#     data_dir = "SimMlTrees304_{}/".format(copy_num)

    print("*********************************************")

    tree_array = merged_file.split("\n")
    
    tree_s = ""
    
#    print(tree_array)
    
#     for i in range(1,305,1):
#         #change the file pattern here
#         file_pat = "/SimMltree_{}.txt".format(i)
#         dir_tree = sim_file_path + file_pat 
#         tree_string = open(dir_tree).read()
#         tree_array.append(tree_string)
#         tree_s += tree_string 

    
#     tree_file = open(sim_file_path+ "/SimMltree_merged.txt", "w+")
#     tree_file.write(tree_s)
#     tree_file.close()    
    
    #CREATE A DATAFRAME TO STORE ALL COMPARISONS OF TREES AND DIFFERENT METHODS OF DISTANCE
    
    col_names = ['Tree_comparison',"Num_First_tree","Num_Second_tree", 'weighted_distance']
    len1 = len(tree_array) - 1
    len2 = int((len1*(len1+1))/2)
#     print(len2)
    tree_data = pd.DataFrame(0, index= range(0,len2), columns = col_names)
    #print(tree_data)
    #print("START!!!")
    
    # Loop to compute pairwise distances between the trees
    count = 0
#     print("count= ",count)
#     print("comparisons num= ",len2)
#     print("dataframe size= ",len(tree_data))
    
    k = 0
    for i in range(k,len(tree_array),1):
        for j in range(k+1,len(tree_array),1):
            #print(str(i+1) + "vs" + str(j+1))
            tree_data.loc[count,'Tree_comparison'] = str(i+1) + "vs" + str(j+1)
            tree_data.loc[count,'Num_First_tree'] = str(i+1)
            tree_data.loc[count,'Num_Second_tree'] = str(j+1)
            
            #Calculate the specified type of tree distance
            if comparison_type == "unweighted_distance":
                tree_data.loc[count,'unweighted_distance'] =unweighted_distance(tree_array[i],tree_array[j])
            elif comparison_type == "weighted_distance":
                
                tree_data.loc[count,'weighted_distance'] = weighted_distance(tree_array[i],tree_array[j])
            else:
                tree_data.loc[count,'euclidean_distance'] = euclidean_distance(tree_array[i],tree_array[j])
                       
            count += 1    
        k += 1;
 
    # Save the results to a CSV file
    tree_data.to_csv('Sim_Trees_distances_comparisons_{}.csv'.format(copy_num), index=False)


    return tree_data
#     #print("END!!")

    
def merge_files(file1,file2,output_path,copy_num):
    t_data = open(file1)
    s1 = str(t_data.read()) + "\n"
    s1 
    t_data2 = open(file2)
    s2 = str(t_data2.read())
    s2 
    s3 = s1 + s2 #Concatenate the contents of the two files

    #writing concatenated contents into a new file
    f= open(output_path + f"{copy_num}.txt","w+")
    f.write(s3)
    f.close()
    return s3




def eric_sum_test(num1,num2):
    s = "you sum for the testing is : " + str(num1 + num2)
    return s
