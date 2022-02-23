# PhanerozoicOceans
[![DOI](https://zenodo.org/badge/455603679.svg)](https://zenodo.org/badge/latestdoi/455603679)
Data and materials from 

> V. Rineau, J, Smyčka, D. Storch, Diversity-dependence is an ubiquitous phenomenon across phanerozoic oceans.

## Installation

This set of scripts has been tested on Ubuntu 20.04 using Rstudio (R version 3.6.3), Python 3.8.10, and openjdk (Java) 17.0.1.

The scripts requires the installation of R and Python 3.

Several packages are also required. For the R scripts, the following lines allow the installation of the required packages.
	
```
install.packages("divDyn")
install.packages("bestNormalize")
install.packages("this.path")
```

The scripts have been tested with the following versions: divDyn 0.8.1, bestNormalize 1.8.2, this.path 0.5.0.

Moreover, the script `4.convergent_cross_maping.R` requires a specific version of the rEDM package (0.7.5) which can be downloaded through devtools ([install devtools on Ubuntu](https://www.digitalocean.com/community/tutorials/how-to-install-r-packages-using-devtools-on-ubuntu-16-04)).

Then the following lines launch the installation of rEDM.

```
library(devtools)
install_github("ha0ye/rEDM")
```

Python dependencies can be easily installed via a terminal using the pip command and the `requirements.txt` file, which contains all the required packages

```
pip install -r requirements.txt
```

Finally, you will need to download and install IDTxl and its dependecies from github: https://github.com/pwollstadt/IDTxl. 
IDTxl requires several packages and Java to work properly: https://github.com/pwollstadt/IDTxl/wiki/Installation-and-Requirements.

## Launch scripts

There are seven scripts. They must be launched in the order of their numbering because each script requires the files produced by the previous one.

1. `1.environmental_datasets.R` - build temperature, carbon, strontium and sulfur time series according to the `Palaeoenvironmental surrogate datasets` SOM part. Result tables are written in the `/datasets/environmental_databases` folder.

Each following script writes results files in `/datasets/taxonomic_databases/taxon_name` for each taxa or `/datasets/taxonomic_databases` for all taxa.

2. `2.building_dataframes.R` - compute diversity, extinction rates, and origination rates according to the `Diversity and turnover rates estimate` SOM part and gather them with the environmental time series in a singe table.
3. `3.cross_correlations.R` - compute cross correlation for each pair of variables and reproduces figures S5 to S13 of the material and methods.
4. `4.convergent_cross_maping.R` - compute convergent cross maping (CCM) analysis for each pair of variables with different lags (figs. S3 and S4). The procedure is given in Material and methods' part entitled `Bivariate causality analyses using Convergent Cross-Mapping`. Result files are written for each taxaset.
5. `5.ccm_networks.py` - reproduces the Figure 1 of the paper following CCM analyses computed in script 4. 
6. `6.Transfer_entropy.py` - computes conditional transfer entropy estimation explained in material and methods, part `Conditional transfer entropy for testing drivers of turnover rates`. The script reproduces figures S145 and S15.
7. `7.var_analyses.py` - compute autoregressive vector modelling (VAR) explained in `Bivariate causality analyses using Convergent Cross-Mapping`. Result files are written for each taxaset.
