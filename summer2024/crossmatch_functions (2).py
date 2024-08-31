import math
import pandas as pd
import numpy as np
import scipy.stats as sc
import networkx as nx

# Function to calculate the binomial coefficient "n choose r"
def choose(n,r):
    f = math.factorial
    return f(n) // f(r) // f(n-r)


#bigN total number of observations
#n total number of trees of one type (in one of the files)

# Function to calculate the cross-match distribution
def crossmatchdist(bigN,n):
    if (bigN%2 == 1):
        return "The number of subjects, bigN, should be even"
    I = int(bigN/2)
    
    dist = pd.DataFrame(index = ["a0", "a1", "a2", "pr"])#creates DataFrame with indices "a0", "a1", "a2", and "pr"
    
    for a1 in range(0,I+1): #loop for possible values of a1
            a2 = (n-a1)/2
            if ((math.floor(a2)==a2)&(a2>=0)):
                    int_a2 = int(a2)
                    a0 = I-(a1+int_a2)
                    
# Notice that a0 = I-(a1+int_a2) = I-(a1 + (n-a1)/2) = I - ( 2a1 + n - a1 )/2 = I - (n+ a1)/2 

                    if (a0>=0): #checking if a0 is non-negative
                            pr = math.factorial(I)/choose(bigN,n) #pr for probability 
                            pr = pr*(2**a1)/(math.factorial(a0)*math.factorial(a1)*math.factorial(int_a2))

                            list1 = [[a0],[a1],[a2],[pr]]

                            dist = pd.concat([dist,pd.DataFrame(list1, index = ["a0", "a1", "a2", "pr"])], axis=1) # Append the values to the DataFrame
                   
                        
    dist.loc["cdf"] = dist.loc["pr"].cumsum() #cdf = cumulative distribution function
    
#    print(dist)
    
    return dist.to_numpy() # Returns the DataFrame as a numpy array

#function to calculate the cross-match distribution and return as DataFrame
def crossmatchdist_df(bigN,n):
    if (bigN%2 == 1):
        return "The number of subjects, bigN, should be even"
    I = int(bigN/2)
    dist = pd.DataFrame(index = [1,2,3,4])

    for a1 in range(0,I+1): #same loop as above for all possible values of a1 
            a2 = (n-a1)/2
            if ((math.floor(a2)==a2)&(a2>=0)):
                    a0 = I-(a1+a2)
                    if (a0>=0):
                            pr = int(math.factorial(I))/choose(bigN,n)
                            pr = int(pr*(2**a1))/(math.factorial(a0)*math.factorial(a1)*math.factorial(a2))
                            list1 = [[a0],[a1],[a2],[pr]]
                            dist = pd.concat([dist,pd.DataFrame(list1)],axis=0)
                            
    
    dist.loc[4] = dist.loc[3].cumsum()
    
    
    return dist

#checking if a matrix is symmetric
def check_symmetric(a, tol=1e-8):
    return np.all(np.abs(a-a.T) < tol)

# Function to perform the cross-match test
def crossmatchtest(z,D):
    if (not check_symmetric(D)): # Check if the distance matrix (D) is symmetric
        raise ValueError("Invalid distance matrix: your distance matrix is not symmetric")
        return NaN
    if (np.any(D<0)): #Check if the distance matrix has any negative values
        raise ValueError("Invalid distance matrix: your distance matrix includes negative values")
        return NaN
    if (len(z) != len(D)): #Check if the vector length matches the matrix size
        raise ValueError("The vector needs to be the same length as the matrix")
        return NaN
        
    
    plainmatrix = 100000*D/D.max() # Normalize the distance matrix
    plainmatrix[np.diag_indices_from(plainmatrix)] = 0 # Set the diagonal to zero
    G = nx.from_numpy_array(D) # Create a graph (G) from the distance matrix
    
    weight = "weight"
    G_edges = G.edges(data=weight, default=1) # Get edges with weights from the graph
    max_weight = 1 + max(w for _, _, w in G_edges) # Calculate the maximum weight.
    InvG = nx.Graph() # Initialize an inverted graph
    edges = ((u, v, max_weight - w) for u, v, w in G_edges) # Invert the edge weights
    InvG.add_weighted_edges_from(edges, weight=weight) # Add weighted edges to the inverted graph
    set1 = nx.max_weight_matching(InvG, maxcardinality=True, weight=weight) #Find the maximum weight matching
        
    list1 = list(set1)  # Convert the matching set to a list
    
    v1 = []
    v2 = []
    
    if(len(list1)>0): #for matched pair
        for i in range(len(list1)):
            first = list1[i][0]
            second = list1[i][1]

            v1.append(first + 1)
            v2.append(second + 1)
            v2.append(first + 1)
            v1.append(second + 1)
    
    min_weight_total = 0

    for h in range(0,len(v1),1): #for summing minimum weights
        if(h%2 == 0):
            min_weight_total = min_weight_total + plainmatrix[v1[h]-1][v2[h]-1]
    
    print("total min weigth is: ", min_weight_total )
    
    v1, v2 = zip(*sorted(zip(v1, v2))) # Sort the pairs
    mt = np.minimum(v1, v2)  # Get the minimum of each pair 
    z0 = z
    
    df = pd.DataFrame({"z0":z0,"mt0":mt}) # Create a DataFrame with z0 and mt0

    df2 = pd.crosstab(index=df['z0'], columns=df['mt0']) # Create a crosstab DataFrame

    sum_total = 0
    for index in df2: #Sum up the values in the crosstab
        if(df2.iloc[0][index] == 1):
            sum_total += 1
    a1 = sum_total
    bigN = len(z0)
    n = sum(z0)
    
    
    if(bigN < 340):        
            dist = crossmatchdist(bigN,n) #Calculate the cross-match distribution
            pval = None
            for j in range(len(dist[1])):
                if(dist[1][j] == a1):
                    pval = dist[3][j] #changed to 3, sum of pr's (??)

    else:
        
        pval = None

    
    m = bigN - n
    Ea1 = (n*m/(bigN-1)) #expected a1 (# of cross matches)
    Va1 = 2*n*(n-1)*m*(m-1)/((bigN-3)*(bigN-1)*(bigN-1)) #variance
    dev = (a1-Ea1)/math.sqrt(Va1) 
    approx = sc.norm.cdf(dev)  #approximate p-value ; should be similar to p-val
    

    final_list = [a1,Ea1,Va1,dev,pval,approx] #Suggestion of change, so it's easier to extract pvalue later
    
    print("--------------------------------------------------------------------------")
    print ("Results of crossmatch test:")  #Suggestion of change for neater output
    print("\n")
    print("a1: ",a1)
    print("Ea1: ",Ea1)
    print("Va1: ",Va1)
    print("dev: ",dev)
    print("pval: ",pval) #reject null if small
    print("approxpval: ",approx)
    print("\n")
    
    
#    final_list = {'a1':a1, 
#                  'Ea1':Ea1,
#                  'Va1':Va1,
#                  'dev':dev,
#                  'pval':pval,
#                  'approxval':approx}   
        
    return final_list