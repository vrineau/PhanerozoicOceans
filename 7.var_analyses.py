#!/usr/bin/env python3
# -*- coding: utf-8 -*-

### V. Rineau, J, Smyčka, D. Storch, Diversity-dependence is an ubiquitous phenomenon across phanerozoic oceans.
### DATA 1.
###
### Data 1: script 5 - 7.var_analyses.py
###
### Computes autoregressive vector modelling (VAR) explained in Material and methods, part
### "Bivariate causality analyses using Convergent Cross-Mapping". 
### Results are written in /datasets/taxonomic_databases.

# Load packages
import numpy as np
import pandas as pd
from statsmodels.tsa.api import VAR
import itertools, sys, os

namefile = ["bivalvia","brachipoda","scleractinia", #list of folders
            "gasteropoda","metazoa","prymnesiophycae",
            "foraminifera","coccolithophoridae","radiolaria"]

namelist = ["Diversity","Extinction rate","Origination rate","Temperature","δ13C","86Sr/87Sr","δ34S"] # List of variables

success = []
successdict = dict()

# For each taxonomic database
for folder in namefile: 
       
    # Change directory
    os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)),
         "datasets","taxonomic_databases",folder))
                           
    env = np.genfromtxt("env_dataframe.csv",delimiter=",", skip_header=1,missing_values="NA",)
        
    # Stating maximum lags
    if folder in namefile[0:5]:
        maxlag = 5
    else:
        maxlag = 10
        
    # Result dataframe
    colnames = []
    for x, y in itertools.permutations([0,1,2,3,4,5,6],2): 
        colnames.append(namelist[x]+" to "+namelist[y])  
    rownames = [*range(maxlag+1)]
    resultsdf = pd.DataFrame(np.nan, index=rownames, columns=colnames)        
    
    # Build table without missing data
    num_rows, num_cols = env.shape
    envsub = env[1:num_rows,[1,2,3,4,5,6,7]]
    envsub = envsub[~np.isnan(envsub).any(axis=1)]
    
    #Launch VAR model
    for x, y in itertools.combinations([0,1,2,3,4,5,6],2): 
                    
        envVAR = pd.DataFrame(envsub[1:num_rows,[x,y]], columns=[namelist[x],namelist[y]])
        
        model = VAR(envVAR)
        results = model.fit(maxlags=maxlag, ic='aic')
        
        original_stdout = sys.stdout # Save a reference to the original standard output
        
        with open('VAR_results.txt', 'a') as f:
            sys.stdout = f # Change the standard output to the created file.
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

    # Save results
    resultsdf.to_csv('VAR_table.csv') #columns: variable links; rows: lags.

