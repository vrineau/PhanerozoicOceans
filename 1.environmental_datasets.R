### V. Rineau, J, Smyčka, D. Storch, Diversity-dependence is an ubiquitous phenomenon across phanerozoic oceans.
### DATA 1.
###
### Data 1: script 1 - 1.environmental_datasets.R
###
### Build temperature, carbon, strontium and sulfur time series according to the 
### material and methods part "Palaeoenvironmental surrogate datasets".  
### Result tables are written in the /datasets/environmental_databases folder.

# Load packages
library(divDyn)
library(this.path)

options(scipen=999) 

# Function for converting ages from one Geological time scale to another
tsconverter <- function(Ao,rosette,i,j) {
  #following Wei and Peleo-Alampay 1993
  column_ages_old <- rosette[i]
  column_ages_old <- column_ages_old[!is.na(column_ages_old)]
  
  To <- column_ages_old[findInterval(Ao,column_ages_old)]
  Bo <- column_ages_old[findInterval(Ao,column_ages_old)+1]
  Tn <- rosette[match(To,rosette[,i]),j]
  Bn <- rosette[match(Bo,rosette[,i]),j]
  An <- (Tn*Ao-Ao*Bn+To*Bn-Tn*Bo)/(To-Bo)
  An
}

# Loess smoothing and bin average (bins defined in stagetable)
binloess <- function(dataset,lissage = 0.09, stagetable) {
  
  mid = dataset$age
  mod <- loess(dataset$mean ~ mid, data = dataset, span = lissage, degree = 2)
  datasetinterpol <- data.frame(mid = stagetable$mid, stg = stagetable$stg)
  dataset2 <- transform(datasetinterpol, interpol = predict(mod, datasetinterpol))
  
  stg <- rep(NA, length(dataset$mean))
  
  for (i in 1:length(dataset$age)) { #nombre de lignes
    for (t in 1:nrow(stagetable)) { #nombre de bins
      nulstg <- FALSE
      if (dataset$age[i]<=stagetable$bottom[t] & dataset$age[i]>=stagetable$top[t]) { # if age is in the bin
        stg[i] <- stagetable$stg[t] #ajout du stg
        nulstg <- TRUE
        break
      }
      if (dataset$age[i]<=stagetable$bottom[nrow(stagetable)] & nulstg == FALSE) { # is lower is minimum
      }
    }
  }
  
  dataset <- cbind(dataset,stg) #ajout de la colonne stg
  dataset <- cbind(dataset,stg.mean=ave(dataset$mean,dataset$stg)) #et les moyennes sans lissage
  
  testdoublon <- (-which(duplicated(dataset[,c(3,4)])))
  
  if (identical(testdoublon, integer(0))==FALSE) {
     dataset <- dataset[-which(duplicated(dataset[,c(3,4)])),] #si doublons, on les supprime
  }
  
  rslt <- na.omit(data.frame(stg = dataset2$stg, mean = dataset2$interpol))
  rslt
  
}

# Load stages and mapping table, and change working directory
data(stages)
setwd(dir = paste(this.dir(), "/datasets/environmental_databases/", sep=""))
rosette <- read.csv("GTS.rosette.csv", sep=";", na.strings="",stringsAsFactors = FALSE)

# Import raw csv files containing environmental time series
C.ogg        <- read.csv("C.ogg.csv", sep=";", na.strings="",stringsAsFactors = FALSE)
T.scotese    <- read.csv("T.scotese.csv", sep=";", na.strings="",stringsAsFactors = FALSE)
S.macarthur  <- read.csv("S.macarthur.csv", sep=";", na.strings="",stringsAsFactors = FALSE)
Sf.prokoph   <- read.csv("Sf.prokoph.csv", sep=";", na.strings="",stringsAsFactors = FALSE)


# Standardisation of column names - extraction and cleaning
C.ogg0 <- na.omit(data.frame(age=as.numeric(as.character(gsub(",",".",C.ogg$age))), 
                             mean=as.numeric(as.character(gsub(",",".",C.ogg$d13C)))))
T.scotese0 <- na.omit(data.frame(age=as.numeric(as.character(gsub(",",".",T.scotese$Age))), 
                                 mean=as.numeric(as.character(gsub(",",".",T.scotese$Mean)))))
S.macarthur0 <- na.omit(data.frame(age=as.numeric(as.character(gsub(",",".",S.macarthur$age))), 
                                   mean=as.numeric(as.character(gsub(",",".",S.macarthur$Mean))))) 
Sf.prokoph0 <- na.omit(data.frame(age=as.numeric(as.character(gsub(",",".",Sf.prokoph$Age..GTS2004.))), 
                                 mean=as.numeric(as.character(gsub(",",".",Sf.prokoph$d34SSSS)))))

# Sort by age
C.ogg0 <- C.ogg0[order(C.ogg0$age),]
T.scotese0 <- T.scotese0[order(T.scotese0$age),]
S.macarthur0 <- S.macarthur0[order(S.macarthur0$age),]
Sf.prokoph0 <- Sf.prokoph0[order(Sf.prokoph0$age),]

# Convert absolute ages to the most recent Geological time Scale 
S.macarthur0$age <- tsconverter(S.macarthur0$age,rosette,3,4)
Sf.prokoph0$age <- tsconverter(Sf.prokoph0$age,rosette,2,4)

S.macarthur0 <- S.macarthur0[!is.na(S.macarthur0$age),]
Sf.prokoph0 <- Sf.prokoph0[!is.na(Sf.prokoph0$age),]

# Smoothing and bin averaging environmental time series following 
# the geological time scale (Paleobiology Database datasets)
C.ogg.dataset        <- binloess(C.ogg0,stages,lissage = 0.01) 
T.scotese.dataset    <- binloess(T.scotese0,stages,lissage = 0.08)
S.macarthur.dataset  <- binloess(S.macarthur0,stages,lissage = 0.05) 
Sf.prokoph.dataset    <- binloess(Sf.prokoph0,stages,lissage = 0.08)

C.ogg.dataset$mean[44:50] <- NA

# Save results for PBDB datasets
write.csv(C.ogg.dataset,"C.ogg.dataset.csv", row.names = FALSE)
write.csv(T.scotese.dataset,"T.scotese.dataset.csv", row.names = FALSE)
write.csv(S.macarthur.dataset,"S.macarthur.dataset.csv", row.names = FALSE)
write.csv(Sf.prokoph.dataset,"Sf.prokoph.dataset.csv", row.names = FALSE)

# Building the stage table for Neptune Database datasets
bins <- 1
xl <- seq(574,0,-bins) #temporal range of the analysis
gstages <- data.frame(sys = NA, system = NA, series = NA, stage = NA, short = NA, bottom = xl, 
                      mid = xl-bins/2, top = xl-bins, dur = bins, stg = 1:length(xl), systemCol = NA, 
                      seriesCol = NA, col = NA)

gstages$top[length(gstages$top)] <- -1

for (i in 1:nrow(gstages)) {
  for (j in 1:nrow(stages)) {
    if (gstages$bottom[i] < stages$bottom[j]) {
      gstages[i,1:5] <- stages[j,1:5]
      gstages[i,11:13] <- stages[j,11:13]
    }
  }
}

# Smoothing and bin averaging environmental time series 
# with 1Myr time bins (Neptune Database datasets)
C.ogg.dataset.micro        <- binloess(C.ogg0,gstages,lissage = 0.01) 
T.scotese.dataset.micro    <- binloess(T.scotese0,gstages,lissage = 0.08)
S.macarthur.dataset.micro  <- binloess(S.macarthur0,gstages,lissage = 0.01) 
Sf.prokoph.dataset.micro    <- binloess(Sf.prokoph0,gstages,lissage = 0.05)

# Save results for Neptune Database datasets
write.csv(C.ogg.dataset.micro,"C.ogg.dataset.micro.csv", row.names = FALSE)
write.csv(T.scotese.dataset.micro,"T.scotese.dataset.micro.csv", row.names = FALSE)
write.csv(S.macarthur.dataset.micro,"S.macarthur.dataset.micro.csv", row.names = FALSE)
write.csv(Sf.prokoph.dataset.micro,"Sf.prokoph.dataset.micro.csv", row.names = FALSE)

# Save mapping table
write.csv(gstages,"gstages.micro.csv", row.names = FALSE)
