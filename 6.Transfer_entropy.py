#!/usr/bin/env python3
# -*- coding: utf-8 -*-

### V. Rineau, J, Smyčka, D. Storch, Diversity-dependence is an ubiquitous phenomenon across phanerozoic oceans.
### DATA 1.
###
### Data 1: script 6 - 6.Transfer_entropy.py
###
### Computes conditional transfer entropy estimation explained in material and methods part entitled 
### "Conditional transfer entropy for testing drivers of turnover rates".
### Results are written in /datasets/taxonomic_databases.

# Load packages
from idtxl.multivariate_te import MultivariateTE
from idtxl.data import Data
import matplotlib.pyplot as plt
import numpy as np
import pandas
import os

# Change directory
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)),
             "datasets","taxonomic_databases"))

namefile = ["bivalvia","brachipoda","scleractinia", #list of folders
            "gasteropoda","metazoa","prymnesiophycae",
            "foraminifera","coccolithophoridae","radiolaria"]

# Load files containing informative time delays
ori_tp = pandas.read_csv("tpvaluesori.csv",delimiter=";")
ext_tp = pandas.read_csv("tpvaluesext.csv",delimiter=";")

successdict = dict()
success = []
result_table = []
i = 0

# Load results object
while i != 8:
    result_table.append(pandas.DataFrame(0, columns = ['Omnibus test', 
                             'Diversity', 'Energy'], index = namefile))
    i += 1

# For each taxonomic dataset
main_folder = os.path.dirname(os.path.abspath(__file__))

for folder in namefile: 
        
    tp_index = -1 # Index of the target variable (0: extinction, 1: origination)
    
    # For extinction and origination as target for CTE analysis
    for tp in [ext_tp,ori_tp]:
        
        tp_index += 1
        taxaindex = tp.index[tp["Unnamed: 0"] == folder][0] #time delay (tp) row number
        
        # For each environmental variable
        for c in [4,5,6,7]:
         
            # Load dataset
            os.chdir(os.path.join(main_folder,folder))
                            
            env = np.genfromtxt("env_dataframe.csv",delimiter=",", 
                                skip_header=1,missing_values="NA",)
        
            # Build new array without missing data
            num_rows, num_cols = env.shape
            envsub = env[1:num_rows,[1,2,3,c]]
            envsub = envsub[~np.isnan(envsub).any(axis=1)]
    
            # Conversion into idtxl format
            data = Data(envsub, dim_order='sp', normalise=False)
    
            # Stating conditioning variables for CTE analysis
            conditionals = [] 
            num_rowssub, num_colssub = envsub.shape
               
            # If informative delays from CCM exists, add them to the 
            # conditioning set for CTE
            for x in [0,3]: #source variable, diversity or environmental variables
                
                if x == 0:
                    y = 1
                else:
                    y = c-3
                    
                if isinstance(tp.iloc[taxaindex,y],str): #if not nan
                    for lag in tp.iloc[taxaindex,y].split():
                        conditionals.append((x,-1*int(lag))) 
                            
            # Initialise analysis object and define settings
            if folder in namefile[0:5]:
                maxlag = 5
            else:
                maxlag = 15
            
            network_analysis = MultivariateTE()
            settings = {'cmi_estimator': "JidtGaussianCMI", #to perform analyses
                        #using the non-linear estimator, change JidtGaussianCMI
                        #by JidtKraskovCMI
                        'max_lag_sources': maxlag,
                        'min_lag_sources': 0, 
                        'verbose': False}
            
            if conditionals:
                settings['add_conditionals'] = conditionals

            # Run analysis
            results = network_analysis.analyse_network(settings=settings, 
                   data=data, targets=[0,1,2]) #targets are div., ext., ori.
                    
            # Fill table with results
            try:
                tabnum = c-5
                
                if tp_index == 0: # Extinction rates
                    results_dict = results.get_single_target(1, fdr=True)
                    
                else: # Origination rates
                    results_dict = results.get_single_target(2, fdr=True)
                    tabnum += 4
                
                # Fill
                result_table[tabnum].iloc[taxaindex,0] = results_dict["omnibus_te"] # Omnibus test results
                
                sources = sorted(list(set([x[0] for x in results_dict["selected_vars_sources"]])))
                
                if 0 in sources:
                    result_table[tabnum].iloc[taxaindex,1] = results_dict["te"][sources.index(0)] # Diversity
                
                if 3 in sources:
                    result_table[tabnum].iloc[taxaindex,2] = results_dict["te"][sources.index(3)] # Environmental variables
                
                result_table[tabnum].iloc[taxaindex,0] = results_dict["omnibus_te"] 
                                
            except:
                pass
            
            if results_dict['omnibus_sign']: #if there is an omnibus value
                for restuple in results_dict['selected_vars_sources']:
                    if restuple[0] == 0: # If diversity signal
                        success.append(folder)                        
                        successdict[folder] = results
                        break

# Build results table
os.chdir(main_folder)
final_table = pandas.concat(result_table, axis=1)
final_table.to_csv("CTE_results.csv")
final_table[final_table < 0] = 0
    
# Build plots    
# Omnibus histogram
Omnibus_table = final_table["Omnibus test"]
Omnibus_table = Omnibus_table.reindex(["prymnesiophycae", "coccolithophoridae", 
                                       "foraminifera", "radiolaria", 
                                       "bivalvia", "gasteropoda", "brachipoda", 
                                       "scleractinia", "metazoa"])
Omnibus_table[Omnibus_table == np.inf] = 2
Omnibus_table.plot.bar(legend=False)

#CTE histogram
fig, axes = plt.subplots(nrows=4, ncols=2)
TE_table = final_table[["Diversity","Energy"]]

i = 0
j = 0

cols = ["Extinction rate","Origination rate"]
rows = ["Temperature","Carbon","Strontium","Sulfur"]

for ax, col in zip(axes[0], cols):
    ax.set_title(col)

for ax, row in zip(axes[:,0], rows):
    ax.set_ylabel(row, rotation=0, size='large')

for p in result_table:
    p[p < 0] = 0
    p = p.reindex(["prymnesiophycae", "coccolithophoridae", 
                                       "foraminifera", "radiolaria", 
                                       "bivalvia", "gasteropoda", "brachipoda", 
                                       "scleractinia", "metazoa"])
    p[["Diversity","Energy"]].plot.bar(ax=axes[j,i],legend=False).set_ylim([0,0.7])
    j += 1
    
    if j == 4:
        j = 0
        i += 1

#Comparative histogram displaying mean of TE diversity-dependence of extinction
#rate and origination rate with all environmental variables as conditioning variables
truc = np.array([result_table[0],result_table[1],result_table[2],result_table[3]])
truc[truc == np.inf] = 2
meantable = truc.mean(axis=0)
stdtable = truc.std(axis=0)
meantablep = pandas.DataFrame(meantable, columns = result_table[0].columns, 
                              index = final_table.index.values)
stdtablep = pandas.DataFrame(stdtable, columns = result_table[0].columns, 
                             index = final_table.index.values)

truc2 = np.array([result_table[4],result_table[5],result_table[6],result_table[7]])
truc2[truc2 == np.inf] = 2
meantable2 = truc2.mean(axis=0)
stdtable2 = truc2.std(axis=0)
meantablep2 = pandas.DataFrame(meantable2, columns = result_table[0].columns, 
                               index = final_table.index.values)
stdtablep2 = pandas.DataFrame(stdtable2, columns = result_table[0].columns, 
                              index = final_table.index.values)

result = pandas.concat([meantablep[["Diversity"]], meantablep2[["Diversity"]]], axis=1)
result = result.set_axis(["extinction rate","origination rate"], axis=1)
result = result.reindex(["prymnesiophycae", "coccolithophoridae", 
                                       "foraminifera", "radiolaria", 
                                       "bivalvia", "gasteropoda", "brachipoda", 
                                       "scleractinia", "metazoa"])

resultstd = pandas.concat([stdtablep[["Diversity"]], stdtablep2[["Diversity"]]], axis=1)
resultstd = resultstd.set_axis(["extinction rate","origination rate"], axis=1)
resultstd = resultstd.reindex(["prymnesiophycae", "coccolithophoridae", 
                                       "foraminifera", "radiolaria", 
                                       "bivalvia", "gasteropoda", "brachipoda", 
                                       "scleractinia", "metazoa"])

result.plot.bar(yerr=resultstd, color = ["red", "green"])
