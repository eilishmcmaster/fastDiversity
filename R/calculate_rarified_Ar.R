calculate_rarified_Ar <- function(gt, site_variable){
  ## calculate rarified allelic richness ##
  # rows are sites and columns are loci 
  count_matrix <- t(sapply(split(seq_along(site_variable), site_variable), 
                           function(idx) colSums(gt[idx, ], na.rm = TRUE))) %>% as.matrix()
  n_matrix <- t(sapply(split(seq_along(site_variable), site_variable), 
                       function(idx) colSums(!is.na(gt[idx, ]))))%>% as.matrix()
  a_matrix <- 2*n_matrix
  #for each locus
  # Determine the smallest n for each locus
  min_n_per_locus <- apply(n_matrix, 2, min)
  n <- min_n_per_locus *2 #minimum allele count per locus
  
  count_matrix2 <- a_matrix-count_matrix
  
  all_ar <- matrix(NA, nrow(n_matrix), ncol(n_matrix)) # matrix to stor all ar values per site and locus
  for(i in 1:ncol(n_matrix)){# for each locus
    for(j in 1:nrow(n_matrix)){# for each site 
      samp <- choose(a_matrix[j,i] - c(count_matrix[j,i], count_matrix2[j,i]), n[i]) /
        choose(a_matrix[j,i], n[i]) # number of alleles at locus i site j - number of reference or alternate alleles observed 
      
      samp[is.na(samp)] <- 0  # Replace any NA values with 0 to handle cases where the probability is undefined
      all_ar[j,i] <- sum(1 - samp)  # Sum the probabilities of sampling each allele to get allelic richness
    }
  }
  return(all_ar)
}