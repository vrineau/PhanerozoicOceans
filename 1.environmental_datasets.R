
library(divDyn) #chargement es package pour la session en cours
library(this.path)

options(scipen=999) 

#FONCTIONS################################################################################################

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

binloess <- function(dataset,lissage = 0.09, stagetable) { #n?cessite que le dataset ait une colonne mean et une colonne age
  
  mid = dataset$age
  mod <- loess(dataset$mean ~ mid, data = dataset, span = lissage, degree = 2)
  datasetinterpol <- data.frame(mid = stagetable$mid, stg = stagetable$stg)
  dataset2 <- transform(datasetinterpol, interpol = predict(mod, datasetinterpol))
  
  stg <- rep(NA, length(dataset$mean))
  
  for (i in 1:length(dataset$age)) { #nombre de lignes
    for (t in 1:nrow(stagetable)) { #nombre de bins
      nulstg <- FALSE
      if (dataset$age[i]<=stagetable$bottom[t] & dataset$age[i]>=stagetable$top[t]) { #si l'age est dans la bin
        stg[i] <- stagetable$stg[t] #ajout du stg
        nulstg <- TRUE
        break
      }
      if (dataset$age[i]<=stagetable$bottom[nrow(stagetable)] & nulstg == FALSE) { #si inf?rieur au plus faible          stg[i] <- stages$stg[nrow(stages)] #ajout du stg
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

#####################################################################################################################
############Traitement des jeux de donn?es: lowess, adaptation aux time bins, cr?ation du env.csv####################
#####################################################################################################################
data(stages)
setwd(dir = paste(dirname(this.dir()), "/datasets/environmental_databases/", sep=""))
rosette <- read.csv("GTS.rosette.csv", sep=";", na.strings="",stringsAsFactors = FALSE)

#import des csv
C.ogg        <- read.csv("C.ogg.csv", sep=";", na.strings="",stringsAsFactors = FALSE)
T.scotese     <- read.csv("T.scotese.csv", sep=";", na.strings="",stringsAsFactors = FALSE)
S.macarthur  <- read.csv("S.macarthur.csv", sep=";", na.strings="",stringsAsFactors = FALSE)
Sf.paytan    <- read.csv("Sf.prokoph.csv", sep=";", na.strings="",stringsAsFactors = FALSE)


#standardisation des noms de colonnes ? extraire et nettoyage
C.ogg0 <- na.omit(data.frame(age=as.numeric(as.character(gsub(",",".",C.ogg$age))), 
                             mean=as.numeric(as.character(gsub(",",".",C.ogg$d13C)))))
T.scotese0 <- na.omit(data.frame(age=as.numeric(as.character(gsub(",",".",T.scotese$Age))), 
                                 mean=as.numeric(as.character(gsub(",",".",T.scotese$Mean)))))
S.macarthur0 <- na.omit(data.frame(age=as.numeric(as.character(gsub(",",".",S.macarthur$age))), 
                                   mean=as.numeric(as.character(gsub(",",".",S.macarthur$Mean))))) 
Sf.paytan0 <- na.omit(data.frame(age=as.numeric(as.character(gsub(",",".",Sf.paytan$Age..GTS2004.))), #GTS 2004->2012 ? convertir
                                 mean=as.numeric(as.character(gsub(",",".",Sf.paytan$d34SSSS)))))

#tri par age
C.ogg0 <- C.ogg0[order(C.ogg0$age),]
T.scotese0 <- T.scotese0[order(T.scotese0$age),]
S.macarthur0 <- S.macarthur0[order(S.macarthur0$age),]
Sf.paytan0 <- Sf.paytan0[order(Sf.paytan0$age),]

#conversion des echelles de temps vers la plus r?cente (2016 ? cette date)
S.macarthur0$age <- tsconverter(S.macarthur0$age,rosette,3,4)
Sf.paytan0$age <- tsconverter(Sf.paytan0$age,rosette,2,4)

S.macarthur0 <- S.macarthur0[!is.na(S.macarthur0$age),]
Sf.paytan0 <- Sf.paytan0[!is.na(Sf.paytan0$age),]

#plot des courbes lissees par LOESS
C.ogg.dataset        <- binloess(C.ogg0,stages,lissage = 0.01) 
T.scotese.dataset    <- binloess(T.scotese0,stages,lissage = 0.08)
S.macarthur.dataset  <- binloess(S.macarthur0,stages,lissage = 0.05) 
Sf.paytan.dataset    <- binloess(Sf.paytan0,stages,lissage = 0.08)

C.ogg.dataset$mean[44:50] <- NA

write.csv(C.ogg.dataset,"C.ogg.dataset.csv", row.names = FALSE)
write.csv(T.scotese.dataset,"T.scotese.dataset.csv", row.names = FALSE)
write.csv(S.macarthur.dataset,"S.macarthur.dataset.csv", row.names = FALSE)
write.csv(Sf.paytan.dataset,"Sf.paytan.dataset.csv", row.names = FALSE)


#####################################################################################################################
############Traitement des jeux de donn?es: lowess, adaptation aux time bins, cr?ation du env.csv####################
#####################################################################################################################

#construction de la stage table
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

#plot des courbes lissees par LOESS - INDICATIF, prise de d?cision dans le bloc suivant
C.ogg.dataset.micro        <- binloess(C.ogg0,gstages,lissage = 0.01) 
T.scotese.dataset.micro    <- binloess(T.scotese0,gstages,lissage = 0.08)
S.macarthur.dataset.micro  <- binloess(S.macarthur0,gstages,lissage = 0.01) 
Sf.paytan.dataset.micro    <- binloess(Sf.paytan0,gstages,lissage = 0.05)

write.csv(C.ogg.dataset.micro,"C.ogg.dataset.micro.csv", row.names = FALSE)
write.csv(T.scotese.dataset.micro,"T.scotese.dataset.micro.csv", row.names = FALSE)
write.csv(S.macarthur.dataset.micro,"S.macarthur.dataset.micro.csv", row.names = FALSE)
write.csv(Sf.paytan.dataset.micro,"Sf.paytan.dataset.micro.csv", row.names = FALSE)

write.csv(gstages,"gstages.micro.csv", row.names = FALSE)
