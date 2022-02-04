#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 11:45:02 2021

@author: valentin
"""

# Import classes
import numpy as np
import pandas as pd
from statsmodels.tsa.api import VAR
import itertools, sys, os

namefile = ["bivalvia","brachipoda","scleractinia", #list of folders
            "gasteropoda","metazoa","prymnesiophycae",
            "foraminifera","coccolithophoridae","radiolaria"]

namelist = ["Diversity","Extinction rate","Origination rate","Temperature","δ13C","86Sr/87Sr","δ34S"]

success = []
successdict = dict()

for folder in namefile: #for each taxonomic database
       
    os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)),
         "datasets","taxonomic_databases",folder))
                           
    env = np.genfromtxt("env_dataframe.csv",delimiter=",", skip_header=1,missing_values="NA",)
        
    #max lags
    if folder in namefile[0:5]:
        maxlag = 5
    else:
        maxlag = 10
        
    #result dataframe
    colnames = []
    for x, y in itertools.permutations([0,1,2,3,4,5,6],2): 
        colnames.append(namelist[x]+" to "+namelist[y])  
    rownames = [*range(maxlag+1)]
    resultsdf = pd.DataFrame(np.nan, index=rownames, columns=colnames)        
    
    #sub-array without missing data
    num_rows, num_cols = env.shape
    envsub = env[1:num_rows,[1,2,3,4,5,6,7]]
    envsub = envsub[~np.isnan(envsub).any(axis=1)]
    
    #VAR model
    for x, y in itertools.combinations([0,1,2,3,4,5,6],2): 
                    
        envVAR = pd.DataFrame(envsub[1:num_rows,[x,y]], columns=[namelist[x],namelist[y]])
        
        model = VAR(envVAR)
        results = model.fit(maxlags=maxlag, ic='aic')
        
        original_stdout = sys.stdout # Save a reference to the original standard output
        
        with open('VAR_results.txt', 'a') as f:
            sys.stdout = f # Change the standard output to the file we created.
            print(results.summary())
            sys.stdout = original_stdout # Reset the standard output to its original value

        for c in [0,1]: #for each column
            
            for r in [*range(1,len(results.params))]: #for each row
            
                if c == 0:
                    if namelist[y] in results.params.index[r]:
                        lagvar = int(results.params.index[r].split(".")[0][1:])
                        relationship = namelist[y] + " to " + namelist[x]
                        resultsdf.iloc[lagvar][relationship] = results.params.iloc[r][c]
            
                if c == 1:
                    if namelist[x] in results.params.index[r]:
                        lagvar = int(results.params.index[r].split(".")[0][1:])
                        relationship = namelist[x] + " to " + namelist[y]
                        resultsdf.iloc[lagvar][relationship] = results.params.iloc[r][c]
    
    resultsdf.loc['mean'] = resultsdf.mean()
    resultsdf.loc['std'] = resultsdf.std()

    resultsdf.to_csv('VAR_table.csv') #columns: variable links; rows: lags.

