import os.path
from AD_Package import panda_plus_symmetry # what I named my own package with the symmetry code
 

# I placed all the file names at the top here so that if they need to be edited then it will
# just be a mass change all throughout the file rather than have to change every single name

# Erics File names
# sim_trees =
# sim_trees_distances = "Sim_Trees_distances_comparisons.csv"

# Aliyahs files
symmetric_matrix = "symmetricMatrix.csv"

# Erics Files
# if not os.path.exists(sim_trees):  # checks if the file already exists
               # if not create the file
# else:  # otherwise it does and print this statement
#     print("The Simulated Data has already been created")

# if not os.path.exists(sim_trees_distances): # checks if the file already exists
                              # if not create the file
# else: # otherwise it does and print this statement
#     print("The Simulated Tree Distances have already been created")


# Aliyahs Files
if not os.path.exists(symmetric_matrix):
    panda_plus_symmetry.panda_plus_symmetry_generate()  # if not create the file
else:
    print("The Symmetric Matrix has already been created")

# Run CrossMatch in python (Eric)


# Save the results to the code
