import math
import pandas as pd
import numpy as np
import scipy.stats as sc
import networkx as nx

def choose(n,r):
    f = math.factorial
    return f(n) // f(r) // f(n-r)


#bigN total number of observations
#n total number of trees of one type (in one of the files)

def crossmatchdist(bigN,n):
    if (bigN%2 == 1):
        return "The number of subjects, bigN, should be even"
    I = int(bigN/2)
    
    dist = pd.DataFrame(index = ["a0", "a1", "a2", "pr"])
    
    for a1 in range(0,I+1):
            a2 = (n-a1)/2
            if ((math.floor(a2)==a2)&(a2>=0)):
                    int_a2 = int(a2)
                    a0 = I-(a1+int_a2)
                    
# Notice that a0 = I-(a1+int_a2) = I-(a1 + (n-a1)/2) = I - ( 2a1 + n - a1 )/2 = I - (n+ a1)/2 

                    if (a0>=0):
                            pr = math.factorial(I)/choose(bigN,n)
                            pr = pr*(2**a1)/(math.factorial(a0)*math.factorial(a1)*math.factorial(int_a2))

                            list1 = [[a0],[a1],[a2],[pr]]

                            dist = pd.concat([dist,pd.DataFrame(list1, index = ["a0", "a1", "a2", "pr"])], axis=1)
                   
                        
    dist.loc["cdf"] = dist.loc["pr"].cumsum() #cdf = cumulative distribution function
    
    print(dist)
    
    return dist.to_numpy()

def crossmatchdist_df(bigN,n):
    if (bigN%2 == 1):
        return "The number of subjects, bigN, should be even"
    I = int(bigN/2)
    dist = pd.DataFrame(index = [1,2,3,4])

    for a1 in range(0,I+1):
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

def check_symmetric(a, tol=1e-8):
    return np.all(np.abs(a-a.T) < tol)

def crossmatchtest(z,D):
    if (not check_symmetric(D)):
        raise ValueError("Invalid distance matrix: your distance matrix is not symmetric")
        return NaN
    if (np.any(D<0)):
        raise ValueError("Invalid distance matrix: your distance matrix includes negative values")
        return NaN
    if (len(z) != len(D)):
        raise ValueError("The vector needs to be the same length as the matrix")
        return NaN
        
    
    plainmatrix = 100000*D/D.max()
    plainmatrix[np.diag_indices_from(plainmatrix)] = 0
    G = nx.from_numpy_array(D)
    
    weight = "weight"
    G_edges = G.edges(data=weight, default=1)
    max_weight = 1 + max(w for _, _, w in G_edges)
    InvG = nx.Graph()
    edges = ((u, v, max_weight - w) for u, v, w in G_edges)
    InvG.add_weighted_edges_from(edges, weight=weight)
    set1 = nx.max_weight_matching(InvG, maxcardinality=True, weight=weight)
        
    list1 = list(set1)
    
    v1 = []
    v2 = []
    
    if(len(list1)>0):
        for i in range(len(list1)):
            first = list1[i][0]
            second = list1[i][1]

            v1.append(first + 1)
            v2.append(second + 1)
            v2.append(first + 1)
            v1.append(second + 1)
    
    min_weight_total = 0

    for h in range(0,len(v1),1):
        if(h%2 == 0):
            min_weight_total = min_weight_total + plainmatrix[v1[h]-1][v2[h]-1]
    
    print("total min weigth is: ", min_weight_total )
    
    v1, v2 = zip(*sorted(zip(v1, v2)))
    mt = np.minimum(v1, v2)  
    z0 = z
    
    df = pd.DataFrame({"z0":z0,"mt0":mt})

    df2 = pd.crosstab(index=df['z0'], columns=df['mt0'])

    sum_total = 0
    for index in df2:
        if(df2.iloc[0][index] == 1):
            sum_total += 1
    a1 = sum_total
    bigN = len(z0)
    n = sum(z0)
    
    
    if(bigN < 340):        
            dist = crossmatchdist(bigN,n)
            pval = None
            for j in range(len(dist[1])):
                if(dist[1][j] == a1):
                    pval = dist[3][j] #changed to 3, sum of pr's (??)

    else:
        
        pval = None

    
    m = bigN - n
    Ea1 = (n*m/(bigN-1))
    Va1 = 2*n*(n-1)*m*(m-1)/((bigN-3)*(bigN-1)*(bigN-1))
    dev = (a1-Ea1)/math.sqrt(Va1)
    approx = sc.norm.cdf(dev)
    

    final_list = [a1,Ea1,Va1,dev,pval,approx] #Suggestion of change, so it's easier to extract pvalue later
    
    print("--------------------------------------------------------------------------")
    print ("Results of crossmatch test:")  #Suggestion of change for neater output
    print("\n")
    print("a1: ",a1)
    print("Ea1: ",Ea1)
    print("Va1: ",Va1)
    print("dev: ",dev)
    print("pval: ",pval)
    print("approxpval: ",approx)
    print("\n")
    
    
#    final_list = {'a1':a1, 
#                  'Ea1':Ea1,
#                  'Va1':Va1,
#                  'dev':dev,
#                  'pval':pval,
#                  'approxval':approx}   
        
    return final_list