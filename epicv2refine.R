#important: beta_matrix is a matrix of beta values (M might also work) with probe names/Illumina IDs in rows and sample names in columns

epicv2refine <- function(beta_matrix) { #whole function takes like 5 mins to run
  #load list of probes to remove, and probes to take a mean of
  load(R-FILE-HERE)
  
  #remove inferior and chromosome 0 probes from beta matrix
  beta2 <- beta_matrix[!(rownames(beta_matrix) %in% PROBES_TO_REMOVE),]
  
  n_replicate_groups <- length(unique(PROBES_TO_MEAN$takemean)) #number of replicate groups = 3329
  
  for (i in 1:n_replicate_groups) {
  current_replicate_group <- PROBES_TO_MEAN$IlmnID[PROBES_TO_MEAN$takemean == (1:n_replicate_groups)[i]] #vector of probes in current replicate group

    #check the beta matrix actually contains all the probes in the group. If it only contains 1 (or 0), no point taking a mean
  temp <- beta2[rownames(beta2) %in% current_replicate_group,]
  if (class(temp)[1] == "matrix") {
    x <- nrow(temp)
    if (x > 0) {
      #print(i) - prints iteration of for loop, useful for debugging but otherwise a bad idea

      #for each sample(/column in beta matrix) replace replicate probe beta values with a mean of the replicate probe beta values. This will result in duplicate rows
      cm <- colMeans(temp, na.rm=TRUE)
      beta2[rownames(beta2) %in% current_replicate_group,] <- matrix(rep(cm, x), nrow = x, byrow = TRUE)
    }
  }
}
  #remove duplicate rows and return result
  beta2 <- unique(beta2)
  return(beta2)
}
