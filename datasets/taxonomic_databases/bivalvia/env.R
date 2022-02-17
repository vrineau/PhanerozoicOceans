library(divDyn)
library(bestNormalize)

#set quorum value for sqs subsampling
qsqs <- 0.8 

#divdyn dataframes
data(stages)
data(stratkeys)
data(keys)
data(tens)

#taxonomic database
taxdb <- read.csv("marine bivalves - pbdb_data.csv", sep = ",", na.strings = "",
                  stringsAsFactors = FALSE, skip = 21) #load taxonomic database
names(taxdb)[6] <- "genus"
taxdb <- taxdb[taxdb$accepted_rank == "genus",]
colnames(tens)[colnames(tens)=="X10"] <- "name"
colnames(stages)[colnames(stages)=="stage"] <- "name"

#database datation
stMin <- categorize(taxdb[ ,"early_interval"], keys$stgInt) #transform stage names in numbers
stMax <- categorize(taxdb[ ,"late_interval"], keys$stgInt) #transform stage names in numbers
stMin <- as.numeric(stMin)
stMax <- as.numeric(stMax)
taxdb$stg <- rep(NA, nrow(taxdb))
stgCondition <- c(
  which(stMax==stMin), # the early and late interval fields indicate the same bin
  which(is.na(stMax))) # or the late_interval field is empty
taxdb$stg[stgCondition] <- stMin[stgCondition] #in these entries, use the bin indicated by the early_interval
taxdbnona <- taxdb[!is.na(taxdb$stg),]

#diversity sqs subsampling and rates calculation
samptax <-binstat(taxdb, tax="genus", bin="stg",coll="collection_no", 
                  ref="reference_no",duplicates=TRUE,indices=TRUE, 
                  noNAStart=TRUE) 

dd <-divDyn(taxdb, bin="stg", tax="genus", noNAStart=TRUE)

sqsquorum <-subsample(taxdbnona, iter=500, q=qsqs, tax="genus", bin="stg", 
                      type="sqs", noNAStart=TRUE)

colnames(sqsquorum)[1] <- "stg"

#environmental datasets
env_path <- paste(substr(getSourceEditorContext()$path,1,35), "datasets/environmental_databases/", sep="")

T.scotese.dataset1    <- read.csv(paste(env_path,"T.scotese.dataset.csv", sep=""),sep = ",", na.strings = "", stringsAsFactors = FALSE) 
C.ogg.dataset1        <- read.csv(paste(env_path,"C.ogg.dataset.csv", sep=""), sep = ",", na.strings = "NA", stringsAsFactors = FALSE) 
S.macarthur.dataset1  <- read.csv(paste(env_path,"S.macarthur.dataset.csv", sep=""), sep = ",", na.strings = "", stringsAsFactors = FALSE) 
Sf.prokoph.dataset1    <- read.csv(paste(env_path,"Sf.prokoph.dataset.csv", sep=""), sep = ",", na.strings = "", stringsAsFactors = FALSE) 

colnames(C.ogg.dataset1)[2] <- "C.ogg"
colnames(T.scotese.dataset1)[2] <- "T.scotese"
colnames(S.macarthur.dataset1)[2] <- "S.macarthur"
colnames(Sf.prokoph.dataset1)[2] <- "Sf.prokoph"

#merging
non_log_env <- Reduce(function(x,y) merge(x,y,by = "stg", all.x = TRUE, all.y = FALSE),
                      list(sqsquorum[,c(1,13,28,29)],
                           T.scotese.dataset1, 
                           C.ogg.dataset1, 
                           S.macarthur.dataset1, 
                           Sf.prokoph.dataset1))

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
env$T.scotese[2:length(env$stg)]    <- diff(env_detrend$T.scotese, differences = 1)
env$S.macarthur[2:length(env$stg)] <- diff(env_detrend$S.macarthur, differences = 1)
env$Sf.prokoph[3:length(env$stg)] <- diff(env_detrend$Sf.prokoph, differences = 2)
env$C.ogg[2:length(env$stg)]    <- diff(env_detrend$C.ogg, differences = 1)

#outlier deletion
env$T.scotese[1] <- NA
env$C.ogg[1] <- NA
env$S.macarthur[1] <- NA
env$Sf.prokoph[1:2] <- NA

#ordernorm transformation and scaling
for (j in 2:ncol(env)) {
  orderNorm_obj <- orderNorm(env[,j])
  env[,j] <- scale(predict(orderNorm_obj))
}

colnames(env) <- c("stage","diversity","extinction","origination","temperature",
                   "carbon","strontium","sulfur")

#save
write.csv(env,"env_dataframe.csv", row.names = FALSE)
save(taxdbnona,env,non_log_env,qsqs,file = "savevar.RData")