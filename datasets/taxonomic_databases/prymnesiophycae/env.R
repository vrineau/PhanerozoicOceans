library(this.path)
library(divDyn)
library(bestNormalize)

#set quorum value for sqs subsampling
qsqs <- 0.9 

#divdyn dataframes
data(stages)
data(stratkeys)
data(keys)
data(tens)

#taxonomic database
taxdb <- read.delim("vrineau_2020-07-30_08-42-55 - diatoms.csv",
                    stringsAsFactors = FALSE)

gstages_path <- paste(this.dir(),"/datasets/environmental_databases/", sep="")

gstages <- read.csv(paste(gstages_path,"gstages.micro.csv", sep=""),
                    sep = ",", na.strings = "", stringsAsFactors = FALSE) 

names(taxdb)[2] <- "genus"
names(taxdb)[19] <- "age"
taxdb <- taxdb[!is.na(taxdb$age),]

#database datation
taxdb$stg <- rep(NA, nrow(taxdb))

#datation par stg bin
for (i in 1:nrow(taxdb)) {
  for (j in 1:nrow(gstages)) {
    if (taxdb$age[i] < gstages$bottom[j] & taxdb$age[i] > gstages$top[j]) {
      taxdb$stg[i] <- gstages$stg[j]
      break
    }
  }
}

taxdb$stg[taxdb$age == 0] <- gstages$stg[length(gstages$stg)]

#diversity sqs subsampling and rates calculation
samptax <-binstat(taxdb, tax="genus", bin="stg",coll="Site", ref="Source", 
                  duplicates=TRUE,indices=TRUE, noNAStart=TRUE) 

dd <-divDyn(taxdb, bin="stg", tax="genus", noNAStart=TRUE) 

sqsquorum <-subsample(taxdb, iter=500, q=qsqs, tax="genus", bin="stg", 
                      type="sqs", noNAStart=TRUE)

colnames(sqsquorum)[1] <- "stg"

#environmental datasets
env_path <- paste(this.dir(),"/datasets/environmental_databases/" sep="")

T.scotese.dataset1    <- read.csv(paste(env_path,"T.scotese.dataset.micro.csv", sep=""),sep = ",", na.strings = "", stringsAsFactors = FALSE)
C.ogg.dataset1        <- read.csv(paste(env_path,"C.ogg.dataset.micro.csv", sep=""), sep = ",", na.strings = "", stringsAsFactors = FALSE)
S.macarthur.dataset1  <- read.csv(paste(env_path,"S.macarthur.dataset.micro.csv", sep=""), sep = ",", na.strings = "", stringsAsFactors = FALSE)
Sf.prokoph.dataset1    <- read.csv(paste(env_path,"Sf.prokoph.dataset.micro.csv", sep=""), sep = ",", na.strings = "", stringsAsFactors = FALSE)

colnames(C.ogg.dataset1)[2] <- "C.veizer"
colnames(T.scotese.dataset1)[2] <- "T.scotese"
colnames(S.macarthur.dataset1)[2] <- "S.macarthur"
colnames(Sf.prokoph.dataset1)[2] <- "Sf.prokoph"

#merging
non_log_env <- Reduce(function(x,y) merge(x,y,by = "stg", all.x = TRUE, all.y = FALSE), 
                      list(sqsquorum[,c(1,13,28,29)],
                           T.scotese.dataset1, 
                           C.ogg.dataset1, 
                           S.macarthur.dataset1, 
                           Sf.prokoph.dataset1
                      ))

#NA deletion
widestts <- c()
for (i in 1:nrow(non_log_env)) {
  tryts <- c()
  for (j in i:nrow(non_log_env)) {
    if (!is.na(non_log_env$divCSIB[j])) {
      tryts <- c(tryts,j) 
    } else {
      break}
  }
  if (length(tryts) > length(widestts)) {widestts <- tryts}
}

non_log_env <- non_log_env[widestts,]

#log
env_detrend <- non_log_env 
env_detrend$divCSIB <- log(non_log_env$divCSIB) 
env <- env_detrend

#remove trend
env$T.scotese    <- c(NA,NA,diff(env_detrend$T.scotese, differences = 2))
env$S.macarthur <- c(NA,NA,diff(env_detrend$S.macarthur, differences = 2))
env$C.veizer    <- c(NA,diff(env_detrend$C.veizer, differences = 1))
env$Sf.prokoph   <- c(NA,NA,diff(env_detrend$Sf.prokoph, differences = 2))

#outlier deletion
env$divCSIB[1] <- NA

#ordernorm transformation and scaling
for (j in 2:ncol(env)) {
  orderNorm_obj <- orderNorm(env[,j])
  env[,j] <- scale(predict(orderNorm_obj))
}

colnames(env) <- c("stage","diversity","extinction","origination","temperature",
                   "carbon","strontium","sulfur")

#save
write.csv(env,"env_dataframe.csv", row.names = FALSE)
taxdbnona <- taxdb
save(taxdbnona,env,non_log_env,qsqs,file = "savevar.RData")

