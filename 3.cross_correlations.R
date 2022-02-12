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
  
  taxdb_path <- paste(dirname(this.dir()), 
                      "datasets/taxonomic_databases/",scr,"/", sep="")
  
  setwd(dir = taxdb_path) #set path
  
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


