import math
import os
import platform

import numpy as np
import pandas as pd


def is_even(val):
    if val % 2 == 0:
        return True
    return False


# takes in the number of elements you want to generate and creates and array of random values
def generate_decimal_array(num_Elements, setState):
    x = []
    if(setState == 0): # don't set a set state
        for i in range(num_Elements):
            array = np.random.uniform()  # https://pynative.com/python-get-random-float-numbers/
            # print(array)
            x.append(array)
    elif(setState >= 1): # set a set state of some value
        np.random.seed(setState)
        for i in range(num_Elements):
            array = np.random.uniform()  # https://pynative.com/python-get-random-float-numbers/
            # print(array)
            x.append(array)

    return x


# Settings that allow all rows and column to be displayed by a dataframe
def see_all(boolean,dataframe):
    x = boolean
    if x:
        pd.set_option('display.max_rows', None)
        pd.set_option('display.max_columns', None)
        pd.set_option('display.width', None)
        pd.set_option('display.max_colwidth', None)
        print(dataframe)



# Code for the actual calculations and generation

#output_folder(string) - file name where the matrix will be saved
# boolean(bool) - will the dataframe be printed
# rep(int) - how many times will the function run, making additional matrix's
# numTrees(int) - how many trees will be in the matrix
# setState(int) - will the same values for the matrix be generated everytime

def dec_matrix_gen(output_folder, boolean, rep, numTrees, setState):
    system = platform.platform() # what platform are you on
    #print(system)

    test_folder_name = output_folder  # folder

    val1 = numTrees ** 2  # the number of trees you want will be squared since I have to make the matrix and it
    # must have a size that will fit all the values generated from the gen decimal array function
    # so each part of the matrix will have sqrt(numtrees **2) elements and there will be sqrt(numtrees **2) groups
    # perfect squares (4, 9, 16, 25, 36, 49, 64, 81, 100, 121, 144, 169, 196, 225, 256, 289, 324, 361, 400, 441, 484, 529, 576, 625, 676, 729, 784, 841, 900)
    # (1..30)^2
    repeat = rep  # how many times you want the data to be generated

    if not os.path.exists(test_folder_name):  # if the folder doesn't exist then make it
        os.makedirs(test_folder_name)

    for iP in range(1, repeat + 1):


        #print("iteration: ", iP, "\n")


        dataframeSetup = []  # array for the index and column names

        for p in range(int(math.sqrt(val1))): #create row and column names for the dataframe
            dataframeSetup.append("Tree" + str(p + 1))  # add the word tree# +1 to the data frame set up (start a tree1)

        ar = generate_decimal_array(val1,setState)
        # print(ar)
        matrix = np.reshape(ar, (int(math.sqrt(val1)), int(math.sqrt(val1))))  # turns the array into a matrix
        # which has sqrt(val) elements broken into sqrt(val) groups


        weightDataFrame = pd.DataFrame(data=np.triu(matrix, k=1), index=dataframeSetup, columns=dataframeSetup) # creates an upper triangle matrix with zeros along a diagonal taking in the matrix

        # Flip the Triangle
        # go through the whole dataframe and copy all the data from one side (1,0) to the other side (0,1) then saves it to
        # the fileName
        for i in range(int(math.sqrt(val1))):
            for a in range(int(math.sqrt(val1))):
                Tree01 = "Tree" + str(a + 1)  # construct the strings form the row name
                Tree02 = "Tree" + str(i + 1)  # construct the strings for the column name

                weightDataFrame.at[Tree01, Tree02] = weightDataFrame.at[Tree02, Tree01]  # insert the data

        # print(weightDataFrame)  # print out the dataframe
        see_all(boolean, weightDataFrame)  # if you want all the information of the data frame to be on

    #system type checker


        if system.__contains__("Windows"):
            fileName = "{folder}\DecimalMatrixGen{num}.csv".format(folder=test_folder_name, num=iP)

        else:  # you are on Mac I hope
            fileName = "{folder}/DecimalMatrixGen{num}.csv".format(folder = test_folder_name, num=iP)

        # print("file_name: ", fileName)
        weightDataFrame.to_csv(fileName) # save data to a csv

#print("runnin")