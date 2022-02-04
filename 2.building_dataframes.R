library(rstudioapi)

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
  
  taxdb_path <- paste(substr(getSourceEditorContext()$path,1,35), 
                      "datasets/taxonomic_databases/",scr,"/", sep="")
  
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

