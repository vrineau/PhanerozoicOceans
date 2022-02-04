library(rstudioapi)
library(rEDM)

scriptlist <- c("bivalvia",
                "brachipoda",
                "scleractinia",
                "gasteropoda",
                "metazoa",
                "prymnesiophycae",
                "foraminifera",
                "radiolaria",
                "coccolithophoridae")

tpvaluesori <- data.frame(diversity=rep(NA, NROW(scriptlist)), 
                       temperature=rep(NA, NROW(scriptlist)), 
                       carbon=rep(NA, NROW(scriptlist)), 
                       strontium=rep(NA, NROW(scriptlist)), 
                       sulfur=rep(NA, NROW(scriptlist)))
rownames(tpvaluesori) <- scriptlist
tpvaluesext <- tpvaluesori
rownames(tpvaluesext) <- scriptlist

for (scr in scriptlist) {
  
  taxdb_path <- paste(substr(getSourceEditorContext()$path,1,35), 
                      "datasets/taxonomic_databases/",scr,"/", sep="")
  
  setwd(dir = taxdb_path) #set path
  
  load("savevar.RData")

  envsub <- env[,2:NCOL(env)]
  envsub <- rawdataset <- as.data.frame(scale(envsub))

  causalvars <- names(envsub)
  effectvars <- causalvars[1:3]

  #simplex : best embedding dimension for ccm analyses #########################
  
  simplex_out <- lapply(causalvars, function(var) {
    simplex(envsub[,var], E = 0:10)
  })
  names(simplex_out) <- causalvars
  
  par(mfrow = c(3, 3))
  for (var in names(simplex_out)) {
    plot(simplex_out[[var]]$E, simplex_out[[var]]$rho, type = "l", 
         xlab = "Embedding Dimension (E)", 
         ylab = "Forecast Skill (rho)", main = var)
  }
  
  best_E <- sapply(simplex_out, function(df) {
    df$E[which.max(df$rho)]
  })

  #convergent cross-maping analyses#############################################
  
  ccm_rho_matrix <- matrix(NA, nrow = length(causalvars), ncol = length(effectvars), dimnames = list(causalvars, effectvars))
  ccm_rho_pvalues_matrix <- matrix(NA, nrow = length(causalvars), ncol = length(effectvars), dimnames = list(causalvars, effectvars))
  ccm_rho_converg <- matrix(NA, nrow = length(causalvars), ncol = length(effectvars), dimnames = list(causalvars, effectvars))
  ccm_tp_matrix <- matrix(NA, nrow = length(causalvars), ncol = length(effectvars), dimnames = list(causalvars, effectvars))
  
  comb <- expand.grid(causalvars, effectvars, stringsAsFactors = FALSE)
  comb <- comb[comb$Var1 != comb$Var2,]
  
  
  par(mfrow = c(3, 3))

  if (scr %in% c("prymnesiophycae","foraminifera","radiolaria",
               "coccolithophoridae")) {tp = -15:0} else {tp = -5:0} #optimal delay (tp) values
  
  #test with several delay values
  for (i in 1:NROW(comb)) {
    vars <- c(comb[i,1], comb[i,2])
    params <- expand.grid(lib_column = vars, target_column = vars, tp = tp) # generate all combinations of lib_column, target_column, tp
    params <- params[params$lib_column != params$target_column, ] # throw out cases where lib == target
    params$E <- best_E[as.vector(params$lib_column)] #embedding dimension
    
    output <- do.call(rbind, lapply(seq_len(NROW(params)), function(i) {
      ccm(envsub, E = params$E[i], lib_sizes = NROW(envsub), 
          random_libs = FALSE, lib_column = params$lib_column[i], 
          target_column = params$target_column[i], 
          tp = params$tp[i], silent = TRUE)}))
    
    plot(output$tp[output$lib_column == comb[i,1]], 
         output$rho[output$lib_column == comb[i,1]], 
         xlab = "tp", ylab = "Cross Map Skill (rho)", col = "green", 
         ylim = c(0,1), type = "l", lwd = 2, main = paste(comb[i,1],"<->",comb[i,2]))
    lines(output$tp[output$lib_column == comb[i,2]], 
          output$rho[output$lib_column == comb[i,2]], col = "purple", 
          lwd = 2)
    legend(x = "topleft", col = c("green", "purple"), lwd = 2, 
           legend = c(paste(comb[i,1]," xmap ",comb[i,2]), paste(comb[i,2],
           " xmap ",comb[i,1])), 
           inset = 0.02, bty = "n", cex = 0.8)
    
    #fill matrices
    ccm_rho_matrix[comb[i,1], comb[i,2]] <- max(output$rho[output$lib_column == comb[i,1]], 
                                                na.rm=TRUE)
    ccm_tp_matrix[comb[i,1], comb[i,2]] <- output$tp[output$rho == ccm_rho_matrix[comb[i,1], 
                                                comb[i,2]] & !is.na(output$rho)]
    
    tp <- output$tp[output$lib_column == comb[i,1]]
    rho <- output$rho[output$lib_column == comb[i,1]]
    
    if (comb[i,2] == "extinction" & comb[i,1] != "origination") {
        tpvaluesext[match(scr, scriptlist), match(comb[i,1], 
                 names(tpvaluesext))] <- paste(tp[rho > 0],collapse=" ")
    } else if (comb[i,2] == "origination" & comb[i,1] != "extinction") {
      tpvaluesori[match(scr, scriptlist), match(comb[i,1], 
                 names(tpvaluesori))] <- paste(tp[rho > 0],collapse=" ")
    }
    
  }
  
  #ccm analyses and significance test with surrogate data
  for (i in 1:NROW(comb)) {
  
    out_temp <- ccm(envsub, E = best_E[comb[i,1]], lib_column = comb[i,1], #ccm
                      target_column = comb[i,2], lib_sizes = NROW(envsub), 
                      tp = ccm_tp_matrix[comb[i,1], comb[i,2]], 
                      replace = FALSE, silent = TRUE)
      
    ccm_rho_matrix[comb[i,1], comb[i,2]] <- out_temp$rho #fill matrix
      
    #significance test
    num_surr <- 1000 #number of runs
    partialdat <- data.frame(envsub[comb[i,1]],envsub[comb[i,2]])
    names(partialdat) <-comb[i,]
    partialdat <- partialdat[complete.cases(partialdat), ]
      
    surr_a <- make_surrogate_data(partialdat[comb[i,2]], method = "ebisuzaki", 
                                    num_surr = num_surr) 
      
    ccm_rho_surr <- numeric(num_surr)
      
    for (j in 1:num_surr) {
      ccm_rho_surr[j] <- ccm(cbind(partialdat[comb[i,1]], surr_a[,j]), 
                               E = best_E[comb[i,1]], lib_column = 1, 
                               tp = ccm_tp_matrix[comb[i,1], comb[i,2]], 
                               target_column = 2, lib_sizes = NROW(envsub), 
                               replace = FALSE, silent = TRUE)$rho
      }
      
    #p-values computation
    ccm_rho_pvalues_matrix[comb[i,1], comb[i,2]] <- (sum(ccm_rho_matrix[comb[i,1], 
                    comb[i,2]] < ccm_rho_surr) + 1) / (length(ccm_rho_surr) + 1)
      
  }
  
  #convergence test
  for (i in 1:NROW(comb)) {
    #cross-maps 
    inv_xmap_no <- ccm(envsub, lib_column = comb[i,1], target_column = comb[i,2], #lib_column cause
                       E = best_E[comb[i,1]], tp = ccm_tp_matrix[comb[i,1], comb[i,2]], silent = TRUE)                       #target_column consequence
    
    no_xmap_inv <- ccm(envsub, lib_column = comb[i,2], target_column = comb[i,1], 
                       E = best_E[comb[i,2]], tp = ccm_tp_matrix[comb[i,1], comb[i,2]], silent = TRUE)
    
    inv_xmap_no_means <- ccm_means(inv_xmap_no) #means for different library sizes
    no_xmap_inv_means <- ccm_means(no_xmap_inv)#moyennes pour différentes library sizes
    
    #convergence test
    kendallconv <- cor.test(inv_xmap_no_means$lib_size, pmax(0, inv_xmap_no_means$rho), method = "kendall")
    if (is.na(kendallconv$estimate) | kendallconv$estimate < 0) {ccm_rho_converg[comb[i,1], comb[i,2]] <- 0.99
    } else  {ccm_rho_converg[comb[i,1], comb[i,2]] <- kendallconv$p.value}
    
  }
  
  ccm_rho_ptot <- ccm_rho_pvalues_matrix
  
  for (i in NROW(ccm_rho_ptot)) {for (j in NCOL(ccm_rho_ptot)) {
    ccm_rho_ptot[i,j] <- max(ccm_rho_pvalues_matrix[i,j],ccm_rho_converg[i,j])
  }}

  #save matrices
  write.table(ccm_rho_matrix,paste("rho_matrix.csv", sep = ""), 
              row.names = TRUE, quote=FALSE,sep=" ")
  write.table(ccm_rho_pvalues_matrix,paste("pval_matrix.csv", sep = ""), 
              row.names = TRUE, quote=FALSE,sep=" ")
  write.table(ccm_rho_converg,paste("converg_matrix.csv", sep = ""), 
              row.names = TRUE, quote=FALSE,sep=" ")
  write.table(ccm_tp_matrix,paste("tp_matrix.csv", sep = ""), 
              row.names = TRUE, quote=FALSE,sep=" ")
  
}

#change directory and save all tp values for TE analysis
setwd("..")
write.csv2(tpvaluesori,file="tpvaluesori.csv")
write.csv2(tpvaluesext,file="tpvaluesext.csv")

