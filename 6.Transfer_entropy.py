#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 11:45:02 2021

@author: valentin
"""

from idtxl.multivariate_te import MultivariateTE
from idtxl.data import Data
import matplotlib.pyplot as plt
import numpy as np
import pandas
import os

os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)),
             "datasets","taxonomic_databases"))

namefile = ["bivalvia","brachipoda","scleractinia", #list of folders
            "gasteropoda","metazoa","prymnesiophycae",
            "foraminifera","coccolithophoridae","radiolaria"]

ori_tp = pandas.read_csv("tpvaluesori.csv",delimiter=";")
ext_tp = pandas.read_csv("tpvaluesext.csv",delimiter=";")

successdict = dict()
success = []
result_table = []
i = 0

#load results object
while i != 8:
    result_table.append(pandas.DataFrame(0, columns = ['Omnibus test', 
                             'Diversity', 'Energy'], index = namefile))
    i += 1

#for each taxonomic dataset
for folder in namefile: 
        
    tp_index = -1 #index of the target variable (0: extinction, 1: origination)
    
    #for extinction and origination as target for CTE analysis
    for tp in [ext_tp,ori_tp]:
        
        tp_index += 1
        taxaindex = tp.index[tp["Unnamed: 0"] == folder][0] #tp row number
        
        #for each environmental variable
        for c in [4,5,6,7]:
         
            #load dataset
            os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)),
             "datasets","taxonomic_databases",folder))
                            
            env = np.genfromtxt("env_dataframe.csv",delimiter=",", 
                                skip_header=1,missing_values="NA",)
        
            #sub-array without missing data
            num_rows, num_cols = env.shape
            envsub = env[1:num_rows,[1,2,3,c]]
            envsub = envsub[~np.isnan(envsub).any(axis=1)]
    
            #conversion in idtxl data object
            data = Data(envsub, dim_order='sp', normalise=False)
    
            conditionals = [] #conditioning variables for CTE analysis
            num_rowssub, num_colssub = envsub.shape
    
                                        
            #if informative delays from CCM exists, add them to the 
            #conditioning set for CTE
            for x in [0,3]: #source variable, diversity or environmental var.
                
                if x == 0:
                    y = 1
                else:
                    y = c-3
                    
                if isinstance(tp.iloc[taxaindex,y],str): #if not nan
                    for lag in tp.iloc[taxaindex,y].split():
                        conditionals.append((x,-1*int(lag))) 
                            
            #initialise analysis object and define settings
            if folder in namefile[0:5]:
                maxlag = 5
            else:
                maxlag = 15
            
            network_analysis = MultivariateTE()
            settings = {'cmi_estimator': "JidtGaussianCMI", #to perform alayses
                        #using the non-linear estimator, change JidtGaussianCMI
                        #by JidtKraskovCMI
                        'max_lag_sources': maxlag,
                        'min_lag_sources': 0, 
                        'verbose': False}
            
            if conditionals:
                settings['add_conditionals'] = conditionals

            #run analysis
            results = network_analysis.analyse_network(settings=settings, 
                   data=data, targets=[0,1,2]) #targets are div., ext., ori.
                    
            #fill results table
            try:
                tabnum = c-5
                
                if tp_index == 0: #ext
                    results_dict = results.get_single_target(1, fdr=True)
                    
                else: #ori
                    results_dict = results.get_single_target(2, fdr=True)
                    tabnum += 4
                
                #remplissage
                result_table[tabnum].iloc[taxaindex,0] = results_dict["omnibus_te"] #omnibus
                
                sources = sorted(list(set([x[0] for x in results_dict["selected_vars_sources"]])))
                
                if 0 in sources:
                    result_table[tabnum].iloc[taxaindex,1] = results_dict["te"][sources.index(0)] #diversity
                
                if 3 in sources:
                    result_table[tabnum].iloc[taxaindex,2] = results_dict["te"][sources.index(3)] #energy
                
                result_table[tabnum].iloc[taxaindex,0] = results_dict["omnibus_te"] #energy
                                
            except:
                pass
            
            if results_dict['omnibus_sign']: #if omnibus TE
                for restuple in results_dict['selected_vars_sources']:
                    if restuple[0] == 0: #if diversity signal
                        success.append(folder)                        
                        successdict[folder] = results
                        break
            
final_table = pandas.concat(result_table, axis=1)
final_table.to_csv("CTE_results.csv")

final_table[final_table < 0] = 0
        
#Omnibus histogram
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
