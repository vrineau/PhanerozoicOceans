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

for (scr in scriptlist) {
  
  #dir(find_root(has_file("2.building_dataframes.R")))
  
  taxdb_path <- paste(dirname(this.dir()), 
                      "/datasets/taxonomic_databases/",scr,"/", sep="")
  
  setwd(dir = taxdb_path) #set path
  
  ll <- parse(file = "env.R")

  #R script loading
  for (i in seq_along(ll)) {
    tryCatch(eval(ll[[i]]), 
      error = function(e) message("catched error", as.character(e)))
    
  }
  
  #cleaning
  rm(list= ls()[!(ls() %in% c("scriptlist","scr"))])
}

