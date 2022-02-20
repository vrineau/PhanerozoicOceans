# PhanerozoicOceans
Data and materials from 

> V. Rineau, J, Smyčka, D. Storch, Diversity-dependence is an ubiquitous phenomenon across phanerozoic oceans.

## Installation


The scripts requires the installation of R and Python 3.

Several packages are also required. For the R scripts, the following lines allow the installation of the required packages.
	
```
install.packages("divDyn")
install.packages("bestNormalize")
install.packages("this.path")
```

Moreover, the script `4.convergent_cross_maping.R` requires a specific version of the rEDM package which can be downloaded through devtools ([install devtools on Ubuntu](https://www.digitalocean.com/community/tutorials/how-to-install-r-packages-using-devtools-on-ubuntu-16-04)).

Then the following lines launch the installation of rEDM.

```
library(devtools)
install_github("ha0ye/rEDM")
```

The Python 3 dependencies can be installed very simply in bash using the file `requirements.txt`, which contains all the required packages, and pip

```
pip install -r requirements.txt
```

Finally, you will need to download and install IDTxl on github: https://github.com/pwollstadt/IDTxl ([see here](https://github.com/pwollstadt/IDTxl/wiki/Installation-and-Requirements)).

## Launch scripts

There are seven scripts. They must be launched in the order of their numbering because each script requires the files produced by the previous one.

1. `1.environmental_datasets.R` - build temperature, carbon, strontium and sulfur time series according to the `Palaeoenvironmental surrogate datasets` SOM part. Result tables are written in the `/datasets/environmental_databases` folder.

Each following script writtes results files in folders held in `/datasets/taxonomic_databases` for each taxa plus the csv files `tpvaluesori.csv` and `tpvaluesext.csv` which compile the informative lag values identified in CCM analyses (script 4).

2. `2.building_dataframes.R` - compute diversity, extinction rate, and origination rates time series according to the `Diversity and turnover rates estimate` SOM part and gather them with the environmental time series in a singe table.
3. `3.cross_correlations.R` - compute cross correlation for each pair of variable.
4. `4.convergent_cross_maping.R` - compute convergent cross maping (CCM) analysis for each pair of variable with different lags. The procedure is given in SOM part entitled `Bivariate causality analyses using Convergent Cross-Mapping`.
5. `5.ccm_networks.py` - reproduces the Figure 1 of the paper following CCM analyses computed in script 4. 
6. `6.Transfer_entropy.py` - computes conditional transfer entropy estimation explained in SOM part entitled `Conditional transfer entropy for testing drivers of turnover rates`.
7. `7.var_analyses.py` - compute autoregressive vector modelling (VAR) explained in `Bivariate causality analyses using Convergent Cross-Mapping`.
