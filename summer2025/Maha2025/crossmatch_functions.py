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

def crossmatchdist_df(bigN, n):
    if (bigN % 2 == 1):
        return "The number of subjects, bigN, should be even"
    I = bigN // 2
    
    # Create empty lists to store values
    a0_list, a1_list, a2_list, pr_list = [], [], [], []
    
    for a1 in range(0, I + 1):
        a2 = (n - a1) / 2
        if a2.is_integer() and a2 >= 0:
            a2 = int(a2)
            a0 = I - (a1 + a2)
            if a0 >= 0:
                a0 = int(a0)
                pr = math.factorial(I) / choose(bigN, n)
                pr *= (2 ** a1) / (math.factorial(a0) * math.factorial(a1) * math.factorial(a2))
                a0_list.append(a0)
                a1_list.append(a1)
                a2_list.append(a2)
                pr_list.append(pr)
    
    # Create DataFrame with proper structure
    dist_df = pd.DataFrame({
        'a0': a0_list,
        'a1': a1_list, 
        'a2': a2_list,
        'pr': pr_list
    })
    
    # Add cumulative distribution
    dist_df['cdf'] = dist_df['pr'].cumsum()
    
    return dist_df

def crossmatchdist(bigN, n):
    if bigN % 2 == 1:
        return "The number of subjects, bigN, should be even"
    
    dist_df = crossmatchdist_df(bigN, n)
    print(dist_df)

    return dist_df.to_numpy()

def check_symmetric(a, tol=1e-8):
    return np.all(np.abs(a - a.T) < tol)

def crossmatchtest(z, D):
    if not check_symmetric(D):
        raise ValueError("Invalid distance matrix: your distance matrix is not symmetric")
    if np.any(D < 0):
        raise ValueError("Invalid distance matrix: your distance matrix includes negative values")
    if len(z) != len(D):
        raise ValueError("The vector needs to be the same length as the matrix")
        
    plainmatrix = 100000 * D / D.max()
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
    
    if len(list1) > 0:
        for i in range(len(list1)):
            first = list1[i][0]
            second = list1[i][1]

            v1.append(first + 1)
            v2.append(second + 1)
            v2.append(first + 1)
            v1.append(second + 1)
    
    min_weight_total = 0

    for h in range(0, len(v1), 1):
        if h % 2 == 0:
            min_weight_total = min_weight_total + plainmatrix[v1[h] - 1][v2[h] - 1]
    
    print("total min weight is: ", min_weight_total)
    
    v1, v2 = zip(*sorted(zip(v1, v2)))
    mt = np.minimum(v1, v2)  
    z0 = z
    
    df = pd.DataFrame({"z0": z0, "mt0": mt})
    df2 = pd.crosstab(index=df['z0'], columns=df['mt0'])

    sum_total = 0
    for index in df2:
        if df2.iloc[0][index] == 1:
            sum_total += 1
    a1 = sum_total
    bigN = len(z0)
    n = sum(z0)
    

    pval = None
    if bigN < 340:        
        dist_df = crossmatchdist_df(bigN, n)
        
        # Calculate one-tailed p-value: sum of probabilities where a1 <= observed a1
        pval = dist_df[dist_df['a1'] <= a1]['pr'].sum()
    else:
        pval = None
    
    # Normal approximation for one-tailed test
    m = bigN - n
    Ea1 = (n * m / (bigN - 1))
    Va1 = 2 * n * (n - 1) * m * (m - 1) / ((bigN - 3) * (bigN - 1) * (bigN - 1))
    dev = (a1 - Ea1) / math.sqrt(Va1)
    
    #We want to use a left-tailed test
    approx = sc.norm.cdf(dev)
    # One-tailed left p-value: probability that a1 <= observed value
    #approx = 1 - sc.norm.cdf(dev) # this is for right-tailed test

    final_list = [a1, Ea1, Va1, dev, pval, approx]
    
    
    return final_list