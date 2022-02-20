### V. Rineau, J, Smyčka, D. Storch, Diversity-dependence is an ubiquitous phenomenon across phanerozoic oceans.
### DATA 1.
###
### Data 1: script 3 - 3.cross_correlations.R
###
### 3.cross_correlations.R - compute cross correlation for each pair of variable. Results are written in 
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

# Loop that load each env.R file from each folder placed in /datasets/taxonomic_databases
# and performs pairwise cross correlation between all pairs of time series.

for (scr in scriptlist) {
  
  taxdb_path <- paste(this.dir(), 
                      "datasets/taxonomic_databases/",scr,"/", sep="")
  
  setwd(dir = taxdb_path) 
  
  load("savevar.RData")

  v <- colnames(env)[2:8]
  
  if (scr %in% c("prymnesiophycae", "foraminifera", "radiolaria", 
               "coccolithophoridae")) {
    lag <- 15
  } else {
    lag <- 5
  }

  vpairs <- t(combn(v, 2))

  par(mfrow = c(3, 3))

  for (vnum in 1:15) {
  
    res_ccf <- ccf(c(env[,vpairs[vnum,1]]), c(env[,vpairs[vnum,2]]),
                   na.action = na.exclude, lag.max=lag, 
                   main=paste(vpairs[vnum,1],"~",vpairs[vnum,2]))
  
  } 
}


