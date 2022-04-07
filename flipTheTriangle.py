import pandas as pd
import numpy as np

dis = pd.read_csv("weightedDistances.csv",header=0,index_col=0)  # loads in data frame recognizing the first row as the names and the first colum as the names
print(dis)

# Test
# print (dis.at['Tree1','Tree2'])
# dis.at['Tree1','Tree2'] = 3333
# print (dis.at['Tree1','Tree2'])
# print (dis)

# go through he whole dataframe and copy all the data from one side (1,0) to the other side (0,1)
for i in range(40):
    for a in range(40):
        Tree02 = "Tree" + str(i+1) # construct the strings for the column name
        Tree01 = "Tree" + str(a+1) # construct the strings form the row name

        dis.at[Tree01,Tree02] = dis.at[Tree02,Tree01] # insert the data

print(dis) # print out the dataframe

dis.to_csv("symmetricMatrix.csv") #

