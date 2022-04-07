# GOAL: filter the trees so that the array only contains 1-40 for trees 1-40
import numpy as np
import pandas as pd

# pd.set_option('display.max_columns', None)  # forces panda to display all the columns
# pd.set_option('display.max_rows', None)

distances_unfiltered = pd.read_csv("Sim_Trees_distances_comparisons.csv")  # loads in data frame

filterArray = [] # array for all the tree filters that will be generated

dataframeSetup = []

final_Num_Trees = 50 # VERY IMPORTANT: this will be the value that is shared through the program the number of trees

# above is the 40 x40 rows(index) and columns that will go into the dataframe
def data_frame_setup(numTrees):
    for p in range(numTrees):
        dataframeSetup.append("Tree" + str(p+1)) # add the word tree# +1 to the data frame set up (start a tree1

data_frame_setup(final_Num_Trees)


weightDataFrame = pd.DataFrame(columns=dataframeSetup, index=dataframeSetup)  # set an empty data frame with 40 columns and 40 rows


for i in range(final_Num_Trees):  # create the filters a put them into a filter array
    tree_filter = (distances_unfiltered["Num_First_tree"] == i) & (distances_unfiltered["Num_Second_tree"] <= final_Num_Trees)  # tree filter is collection of boolean values literally a filter filters the trees into boolean values
    filterArray.append(tree_filter)


for i in range(final_Num_Trees):
    zeroArray = [0]  # this stores the zeros that will go in the data frame for comparisons that already were made
    # like Tree 1 Vs Tree 1 or Tree 2 vs Tree 1 (a comparison made in the previous row
    # also, here so that the zeros sorted within can be erased and then refilled depending on the value of i
    if i != (final_Num_Trees-1): # all filters except the last one
        distances_filtered = distances_unfiltered[filterArray[i+1]]  # filter 1 (filter 1 is at index 1 so do i + 1)

        weightedDistancesTemp = distances_filtered["weighted_distance"]  # a series pulling the one column from distances filtered


        for k in range(i):
            zeroArray = np.append(0, zeroArray) # add i-1 zeros to the array

        appendZero = np.append(zeroArray, weightedDistancesTemp.values)  # insert a zeros at the front of the array (to fill it to 40)


        # Place the data for each tree into the weightDataFrame

        treeRow = "Tree" + str(i+1)  # there is no tree zero (Tree0) so start at tree i + 1
        weightDataFrame.loc[treeRow] = appendZero  # insert array into the first row of dataframe

    else:
        for x in range(final_Num_Trees-1):  # fill the zero array with zeros
            zeroArray = np.append(0, zeroArray)


        # Place the data for each tree into the weightDataFrame
        weightDataFrame.loc["Tree" + str(final_Num_Trees)] = zeroArray  # insert zeros into the last row of dataframe (all comparisons would have already been made)


# Flip the Triangle

# go through the whole dataframe and copy all the data from one side (1,0) to the other side (0,1)
for i in range(final_Num_Trees):
    for a in range(final_Num_Trees):
        Tree01 = "Tree" + str(a+1)  # construct the strings form the row name
        Tree02 = "Tree" + str(i + 1)  # construct the strings for the column name

        weightDataFrame.at[Tree01, Tree02] = weightDataFrame.at[Tree02, Tree01]  # insert the data

print(weightDataFrame)  # print out the dataframe

weightDataFrame.to_csv("symmetricMatrix.csv")
