import os.path

import pandas as pd

from AD_Package import panda_plus_symmetry # what I named my own package with the symmetry code
 

# I placed all the file names at the top here so that if they need to be edited then it will
# just be a mass change all throughout the file rather than have to change every single name

# Erics File names and variables
# sim_trees =
# sim_trees_distances = "Sim_Trees_distances_comparisons.csv"
#sim_tree_dis_DF = # this is the variable for the dataframe for the distances

# Aliyahs files names and variables
symmetric_matrix = "symmetricMatrix4.csv"

# Erics Files
# if not os.path.exists(sim_trees):  # checks if the file already exists

               # if not create the file
# else:  # otherwise it does and print this statement
#     print("The Simulated Data has already been created")

# if not os.path.exists(sim_trees_distances): # checks if the file already exists
# sim_tree_dis_Matrix
                              # if not create the file
# else: # otherwise it does and print this statement
#     print("The Simulated Tree Distances have already been created")


#sim_tree_dis_DF = pd.read_csv("Sim-Trees-distances-comparisons-1.csv")
# Aliyahs Files
if not os.path.exists(symmetric_matrix):
    panda_plus_symmetry.panda_plus_symmetry_generate(sim_tree_dis_DF, 50).to_csv(symmetric_matrix) # this won't work until a dataframe is passed into it so comment this line out  basically all Aliyahs code if your working on something else here

else:
    print("The Symmetric Matrix has already been created")

# Run CrossMatch in python (Eric)


# Save the results to the code
