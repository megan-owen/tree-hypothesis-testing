import numpy as np
import pandas as pd
from newest_Eric_crossmatchtest import crossmatchtest
from AD_Package.decimal_matrix_generator import dec_matrix_gen
import platform

system = platform.platform()
#print(system)

#SETTINGS
setState_code = 0 # will the same matrixes be generated
boolean_show = False # will the whole
#results = open("Python_Crossmatch_Results1.csv", "w") # open a file
soop= r"/"

if system.__contains__("Windows"): #change slash if on windows
    soop = r"\""

# in_folder = dir containing matrix's
#file_name - the name of the file in the dir

def cross_test(in_folder, file_name): # runs crossmatch test
    test_folder = in_folder
    df = pd.read_csv("{folder}{slash}{file}.csv".format(folder=test_folder, file=file_name, slash = soop ), header=0,
                     index_col=0)  # import data
    df_ar = pd.DataFrame.to_numpy(df)  # turn dataframe  into a numpy array
    z_data = np.reshape(df_ar, ((int(len(df_ar))), (int(len(df_ar)))))  # takes the data and put them into a matrix
    vec = np.concatenate((np.repeat([0], len(df) / 2), np.repeat([1], len(df) / 2))).tolist()  # creates vector
    return (crossmatchtest(vec, z_data))  # run the crossmatch test


# p - folder number containing the matrix
# rep number of files you want to generate from the function dec_matrix_gen , this is also the num of matrix's in a dir
# numtrees - the number of trees contained within each matrix generated
# NOTE : CURRENTLY the matrix generators are commented out so existing matrices are needed
def printDat (p, rep, numtrees): # prints data to screen
    print(  "-------------------------------------TESTING {val}----------------------------------------------------------------\n".format(val = p))
    #dec_matrix_gen("testFolder{val}".format(val = p), False, rep, numtrees, 0)
    for i in range(1, rep+1): # go through the whole folder rep = the number of files in a folder, +1 to do all files
        print("---------------------------------------------------------------------------------------\n")
        print("crossmatchtest run for DecimalMatrixGen{num}".format(num=i))
        print(cross_test("testFolder{val}".format(val = p), "DecimalMatrixGen{num}".format(num=i)))
        print("---------------------------------------------------------------------------------------\n")

def printDat2 (p, rep, numtrees): # print to a file
    print("calculating crossmatch for files in testFolder{val}... ".format(val= p))
    #dec_matrix_gen("testFolder{val}".format(val = p), False, rep, numtrees, 0)
    data_f = pd.DataFrame(columns=["a1","Ea1","Va1","pval","dev","approxval" ]) # have to reset and empty it for next time
    results = open("Python_Crossmatch_Results"+str(p)+".csv", "w")  # open a file
    #print("Python_Crossmatch_Results"+str(p)+".csv")
    for i in range(1, rep + 1):  # go through the whole folder rep = the number of files in a folder
        diction = cross_test("testFolder{val}".format(val=p), "DecimalMatrixGen{num}".format(num=i))
        # print(diction)
        data_f.loc[i] = diction

    data_f.index = ["DecimalMatrixGen1","DecimalMatrixGen2","DecimalMatrixGen3","DecimalMatrixGen4","DecimalMatrixGen5","DecimalMatrixGen6","DecimalMatrixGen7","DecimalMatrixGen8","DecimalMatrixGen9","DecimalMatrixGen10"]
    data_f.to_csv(results)
    results.close() # close the document
    print("COMPLETED crossmatch for files in testFolder{val} \n".format(val=p))
    #print(data_f)


# Comparing Python with R
def compare_R_Python(rCSV,pyCSV):
    diff_count = 0
    f= pd.read_csv(rCSV)
    #print("f\n",f)
    g = pd.read_csv(pyCSV, index_col= 0)
    #print("g\n",g)


    if not (f.size == g.size): # if the sizes are diffrent
        #print("f_R.size:", f.size)
        #print("g_PY.size:", g.size)
        #print("The dataframes are not the same size")
        return False

    listOcol = ["a1","Ea1","Va1","pval","dev","approxval"]
    # print(f.at["DecimalMatrixGen1","a1"])

    for j in range(6): # last element is approx so len of list is 6
        for i in range(1, len(f)): # run through the dataframe to see if there are any differing values
            row = "DecimalMatrixGen"+ str(i)
            if abs(f.at[row,listOcol[j]] - g.at[row,listOcol[j]]) > 0.000001: # they differ Erics code reformatted
                #print("there is a difference between f: ", f[i], " and ", g[i])
                diff_count += 1
                print("DIFFERENCE FOUND:")
                print("in f_Py at ("+ row +","+ listOcol[j] + ") there is "+ str(f.at[row,listOcol[j]]))
                print("in g_R at (" + row + ", " + listOcol[j] + ") there is " + str(g.at[row, listOcol[j]])+ "\n")
                #print("diff+1")

        #print("checked", listOcol[j])


    if (diff_count > 0): # if theres a diffrence
        #print("Dataframes are not equal\n")
        return False

    #print("the Dataframes are equal")
    return True



#Test Cases

#create mock matrixes 1 and 2 are the same 3 is diffrent
# dec_matrix_gen("pizza1",False,1,4,1)
# dec_matrix_gen("pizza2",False,1,4,1)
# dec_matrix_gen("pizza3",False,1,4,0)

#write the operations to the screen
#printDat(1, 10, 4)
# printDat(2, 10, 10)
# printDat(3, 10, 30)
# printDat(4, 10, 40)
# printDat(5, 10, 80)
# printDat(6, 10, 100)
# printDat(7, 10, 200)
# printDat2(8,10,250)
# printDat2(9,10,304)

#Write the operations to a file
# printDat2(1, 10, 4)
# printDat2(2, 10, 10)
# printDat2(3, 10, 30)
# printDat2(4,10,40)
# printDat2(5,10,80)
# printDat2(6,10,100)
# printDat2(7,10,200)
# printDat2(8,10,250)
# printDat2(9,10,304)


#TESTS fully written out
# FORMAT : dec_matrix_gen(output_folder, boolean, rep, numTrees, setState)

#
# #Test 1 - Generating ten,4 tree decimal matrix
# print("-------------------------------------TESTING 1----------------------------------------------------------------\n")
# dec_matrix_gen("testFolder1", boolean_show,10, 4, setState_code)
# for i in range (1,11):
#     print("---------------------------------------------------------------------------------------\n")
#     print("crossmatchtest run for DecimalMatrixGen{num}".format(num = i))
#     print(cross_test("testFolder1", "DecimalMatrixGen{num}".format(num = i)))
#     print("---------------------------------------------------------------------------------------\n")
#
# #
# # Test 2 - Generating ten, 10 tree decimal matrix
# print("-------------------------------------TESTING 2----------------------------------------------------------------\n")
# dec_matrix_gen("testFolder2", boolean_show, 10,10,setState_code)
# for i in range (1,11):
#     print("---------------------------------------------------------------------------------------\n")
#     print("crossmatchtest run for DecimalMatrixGen{num}".format(num = i))
#     print(cross_test("testFolder2", "DecimalMatrixGen{num}".format(num = i)))
#     print("---------------------------------------------------------------------------------------\n")
#
# # Test 3 - Generating ten, 30 tree decimal matrix
# print("-------------------------------------TESTING 3----------------------------------------------------------------\n")
# dec_matrix_gen("testFolder3", boolean_show, 10,30,setState_code)
# for i in range (1,11):
#     print("---------------------------------------------------------------------------------------\n")
#     print("crossmatchtest run for DecimalMatrixGen{num}".format(num = i))
#     print(cross_test("testFolder3", "DecimalMatrixGen{num}".format(num = i)))
#     print("---------------------------------------------------------------------------------------\n")
#
#
# # Test 4 - Generating ten, 40 tree decimal matrix
# print("-------------------------------------TESTING 4----------------------------------------------------------------\n")
# dec_matrix_gen("testFolder4", boolean_show, 10,40, setState_code)
# for i in range (1,11):
#     print("---------------------------------------------------------------------------------------\n")
#     print("crossmatchtest run for DecimalMatrixGen{num}".format(num = i))
#     print(cross_test("testFolder4", "DecimalMatrixGen{num}".format(num = i)))
#     print("---------------------------------------------------------------------------------------\n")
#
# # Test 5 - Generating ten, 80 tree decimal matrix
# print("-------------------------------------TESTING 5----------------------------------------------------------------\n")
# dec_matrix_gen("testFolder5", boolean_show, 10,80,setState_code)
# for i in range (1,11):
#     print("---------------------------------------------------------------------------------------\n")
#     print("crossmatchtest run for DecimalMatrixGen{num}".format(num = i))
#     print(cross_test("testFolder5", "DecimalMatrixGen{num}".format(num = i)))
#     print("---------------------------------------------------------------------------------------\n")
#
#
# # Test 6 - Generating ten, 100 tree decimal matrix
# print("-------------------------------------TESTING 6----------------------------------------------------------------\n")
# dec_matrix_gen("testFolder6", boolean_show, 10,100,setState_code)
# for i in range (1,11):
#     print("---------------------------------------------------------------------------------------\n")
#     print("crossmatchtest run for DecimalMatrixGen{num}".format(num = i))
#     print(cross_test("testFolder6", "DecimalMatrixGen{num}".format(num = i)))
#     print("---------------------------------------------------------------------------------------\n")
#
# # Test 7 - Generating ten, 200 tree decimal matrix
# print("-------------------------------------TESTING 7----------------------------------------------------------------\n")
# dec_matrix_gen("testFolder7", boolean_show, 10,200,setState_code)
# for i in range (1,11):
#     print("---------------------------------------------------------------------------------------\n")
#     print("crossmatchtest run for DecimalMatrixGen{num}".format(num = i))
#     print(cross_test("testFolder7", "DecimalMatrixGen{num}".format(num = i)))
#     print("---------------------------------------------------------------------------------------\n")
#
#
# # Test 8 - Generating ten, 250 tree decimal matrix
# print("-------------------------------------TESTING 8----------------------------------------------------------------\n")
# dec_matrix_gen("testFolder8", boolean_show, 10,250,setState_code)
# for i in range (1,11):
#     print("---------------------------------------------------------------------------------------\n")
#     print("crossmatchtest run for DecimalMatrixGen{num}".format(num = i))
#     print(cross_test("testFolder8", "DecimalMatrixGen{num}".format(num = i)))
#     print("---------------------------------------------------------------------------------------\n")
#
#
# # Test 9 - Generating ten, 304 tree decimal matrix
# print("-------------------------------------TESTING 9----------------------------------------------------------------\n")
# dec_matrix_gen("testFolder9", boolean_show, 10,304,setState_code)
# for i in range (1,11):
#     print("---------------------------------------------------------------------------------------\n")
#     print("crossmatchtest run for DecimalMatrixGen{num}".format(num = i))
#     print(cross_test("testFolder9", "DecimalMatrixGen{num}".format(num = i)))
#     print("---------------------------------------------------------------------------------------\n")
#
# print("----------------------------------COMPLETE-----------------------------------------------------\n")



#RUN THE COMPARISON
# change the 10 for whatever value so that you can run through that number of files to compare
for i in range(1,10): # remember the second val is not included so there will be no Crossmathc_Results10.txt
    file1 = open("R_Crossmatch_Results"+str(i)+".csv")
    file2 = open("Python_Crossmatch_Results"+str(i)+".csv")

    print("comparing " + "R_Crossmatch_Results"+str(i)+".csv to "+"Python_Crossmatch_Results"+str(i)+".csv \n" )
    if(compare_R_Python(file1,file2)):
        print("for R_Crossmatch_Results"+str(i)+".csv vs."+ "Python_Crossmatch_Results"+str(i)+".csv:" + "\nthe Dataframes are equal")
    else:
        print("for R_Crossmatch_Results"+str(i)+".csv vs."+ "Python_Crossmatch_Results"+str(i)+".csv:" + "\nthe Dataframes are NOT equal")
    print("COMPARISON COMPLETE\n")


#Initial Generating Settings Brief
# This Code by default can generate 9 Test Folders
# Each Contains 10 decimal matrix's
# 9 "Result" Files
# Each Conatining 10 rows from a dataframe 6 columns "a1","Ea1","Va1","pval","dev","approxval"