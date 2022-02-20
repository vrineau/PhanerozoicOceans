### V. Rineau, J, Smyčka, D. Storch, Diversity-dependence is an ubiquitous phenomenon across phanerozoic oceans.
### DATA 1.
###
### Data 1: script 2 - 2.building_dataframes.R
###
### Compute diversity, extinction rate, and origination rates time series according to the 
### "Diversity and turnover rates estimate" material and methods part and gather them with 
### the environmental time series in a singe table. Results are written in 
### /datasets/taxonomic_databases.

library(this.path)

scriptlist <- c("bivalvia",
                "brachipoda",
                "scleractinia",
                "gasteropoda",
                "metazoa",
                "prymnesiophycae",
                "foraminifera",
                "radiolaria",
                "coccolithophoridae")

# Loop that run each env.R file from each folder placed in /datasets/taxonomic_databases
# The resulting tables are named 'env_dataframe.csv' and contain the following 
# columns (each column is a detrended time series) :
# * diversity
# * extinction rates
# * origination rates
# * temperature
# * carbon (Delta13C)
# * strontium (86Sr/87Sr)
# * sulfur (Delta34S)

for (scr in scriptlist) {
    
  taxdb_path <- paste(this.dir(), 
                      "/datasets/taxonomic_databases/",scr,"/", sep="")
  
  # Set path
  setwd(dir = taxdb_path)
  
  ll <- parse(file = "env.R")

  # R script loading
  for (i in seq_along(ll)) {
    tryCatch(eval(ll[[i]]), 
      error = function(e) message("catched error", as.character(e)))
    
  }
  
  # Cleaning
  rm(list= ls()[!(ls() %in% c("scriptlist","scr"))])
}

