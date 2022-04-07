# GOAL: filter the trees so that the array only contains 1-40 for trees 1-40
import numpy as np
import pandas as pd

# pd.set_option('display.max_columns', None)  # forces panda to display all the columns
# pd.set_option('display.max_rows', None)

distances_unfiltered = pd.read_csv("Sim_Trees_distances_comparisons.csv")  # loads in data frame
#print(distances_unfiltered)

filterArray = [] # array for all the tree filters that will be generated
dataframeSetup = ['Tree1','Tree2','Tree3','Tree4','Tree5','Tree6', 'Tree7','Tree8','Tree9','Tree10',
                  'Tree11','Tree12','Tree13','Tree14','Tree15','Tree16','Tree17','Tree18','Tree19','Tree20',
                  'Tree21','Tree22','Tree23','Tree24','Tree25','Tree26','Tree27','Tree28','Tree29','Tree30',
                  'Tree31','Tree32','Tree33','Tree34','Tree35','Tree36','Tree37','Tree38','Tree39','Tree40']
# above is the 40 x40 rows(index) and columns that will go into the dataframe

weightDataFrame = pd.DataFrame(columns=dataframeSetup, index=dataframeSetup)  # set an empty data frame with 40 columns and 40 rows


#  weightDataFrame.Tree1['Tree2'] = '1' # inserts data value 1 into column Tree1 and row Tree2
#print(weightDataFrame)


for i in range(40):  # create the filters a put them into a filter array
    tree_filter = (distances_unfiltered["Num_First_tree"] == i) & (distances_unfiltered["Num_Second_tree"] <= 40)  # tree filter is collection of boolean values literally a filter filters the trees into boolean values
    filterArray.append(tree_filter)

# for filters in filterArray:  # go through all the filters in the filterArray plug them in to weighted Distances
#
#     distances_filtered = distances_unfiltered[filters]  # filter the unfiltered data
#
#     weightedDistancesTemp = (distances_filtered["weighted_distance"])
#
#     print(weightedDistancesTemp)


for i in range(40):
    zeroArray = [0]  # this stores the zeros that will go in the data frame for comparisons that alreadyt were made like Tree 1 Vs Tree 1 or Tree 2 vs Tree 1 (a comparsion made in the previous row
    # also here so that the zeros sotored within can be erased and then refilled depending on the value of i
    if i != 39:# since the whole program is stating at index position 1 for the filters when we get to 40 there will be a call for a 41st filter which doesnt exist so to fix that we run the program up until that point
        # then since all the comparisons who have already been made comparing tree 40 to everything else the last row will be filled with zeros
        distances_filtered = distances_unfiltered[filterArray[i+1]]  # filter 1 (filter 1 is at index 1 so do i + 1)
        #print(distances_filtered)

        weightedDistancesTemp = distances_filtered["weighted_distance"]  # a series

        #print(weightedDistancesTemp.values) # pulls out the valuse of a series and turns it into an array

        for k in range(i):
            zeroArray = np.append(0, zeroArray) # add i-1 zeros to the array

        appendZero = np.append(zeroArray, weightedDistancesTemp.values)  # insert a zeros at the front of the array (to fill it to 40)
        #print(appendZero)

        # Place the data for each tree into the weightDataFrame

        treeRow = "Tree" + str(i+1) #there is no tree zero so start at tree i + 1
        weightDataFrame.loc[treeRow] = appendZero  # insert array into the first row of dataframe

    else:
        for x in range(39): # here we start at x = 0 unlike in the part of the if statement for the filters so we are going to 39
            zeroArray = np.append(0,zeroArray)
        # print(appendZero)

        # Place the data for each tree into the weightDataFrame
        weightDataFrame.loc["Tree40"] = zeroArray  # insert array into the last row of dataframe

print(weightDataFrame)





#Save the weight dataframe to a csv
weightDataFrame.to_csv('weightedDistances.csv')

