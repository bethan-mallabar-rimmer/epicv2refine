#important: beta_matrix is a matrix of beta values (M might also work) with probe names/Illumina IDs in rows and sample names in columns
#by default the function calculates variance of probe signal within replicate groups, and removes outliers with very high variance
#set remove_high_var = FALSE to take a mean of replicates with high variance, instead of removing them

epicv2refine <- function(beta_matrix, remove_high_var = TRUE, print_loop = FALSE) { #whole function takes like 5 mins to run
  #load list of probes to remove, and probes to take a mean of
  load(url('https://github.com/bethan-mallabar-rimmer/epicv2refine/raw/main/FINAL_PROBES.RData'))
  
  if (substr(PROBES_TO_MEAN$takemean[1],1,4) == "mean") {
    #the PROBES_TO_MEAN file I uploaded to github has replicate groups labelled with the format 1, 2, 3... but sometimes it randomly
    #reverts back to a previous format 'mean1', 'mean2', 'mean3' so this bit of code hopefully fixes that
    PROBES_TO_MEAN$takemean <- substr(PROBES_TO_MEAN$takemean, 5, nchar(PROBES_TO_MEAN$takemean))
    PROBES_TO_MEAN$takemean <- as.numeric(PROBES_TO_MEAN$takemean)
  }
  
  #remove inferior and chromosome 0 probes from beta matrix
  beta2 <- beta_matrix[!(rownames(beta_matrix) %in% PROBES_TO_REMOVE),]
  
  n_replicate_groups <- length(unique(PROBES_TO_MEAN$takemean)) #number of replicate groups = 3329
  
  if (remove_high_var == TRUE) {
    #calculate variance; probes with high variance will be removed as it may not be valid
    #to either take a mean or select a random replicate from the group
    var_replicate_groups <- data.frame(meangroup = 1:n_replicate_groups,
                                       var = rep(0, n_replicate_groups),
                                       otlr = rep(FALSE, n_replicate_groups))
  }
  
  
  for (i in 1:n_replicate_groups) {
    current_replicate_group <- PROBES_TO_MEAN$IlmnID[PROBES_TO_MEAN$takemean == (1:n_replicate_groups)[i]] #vector of probes in current replicate group
    
    #check the beta matrix actually contains all the probes in the group. If it only contains 1 (or 0), no point taking a mean
    temp <- beta2[rownames(beta2) %in% current_replicate_group,]
    if (class(temp)[1] == "matrix") {
      x <- nrow(temp)
      if (x > 0) {
        if  (print_loop == TRUE) {
        print(i) #prints iteration of for loop, useful for debugging but otherwise a bad idea
        }
        
        if (remove_high_var == TRUE) {
          #calculate variance of this replicate group (including all replicates in all samples)
          var_replicate_groups$var[i] <- var(as.vector(beta2[rownames(beta2) %in% current_replicate_group,]), na.rm = TRUE)
        }
        
        #for each sample(/column in beta matrix) replace replicate probe beta values with a mean of the replicate probe beta values. This will result in duplicate rows
        cm <- colMeans(temp, na.rm=TRUE)
        beta2[rownames(beta2) %in% current_replicate_group,] <- matrix(rep(cm, x), nrow = x, byrow = TRUE)

      }
    }
  }
  
  if (remove_high_var == TRUE) {
    #mark variance outliers
    upper_bound <- quantile(var_replicate_groups$var)[4] + (1.5 * IQR(var_replicate_groups$var))
    var_replicate_groups$otlr[var_replicate_groups$var >= upper_bound] <- TRUE
    
    #remove outliers with high variance
    #no need to exclude outliers with low variance, as less is better
    for (i in var_replicate_groups$meangroup[var_replicate_groups$otlr == TRUE]) {
      outlier_replicate_group <- PROBES_TO_MEAN$IlmnID[PROBES_TO_MEAN$takemean == i]
      beta2 <- beta2[-(which(rownames(beta2) %in% outlier_replicate_group)),]
    }
  }
  
  #remove duplicate rows and return result
  beta2 <- unique(beta2)
  return(beta2)
}
