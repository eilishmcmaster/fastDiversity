faststats <- function(gt, genetic_group_variable, site_variable, minimum_n=3, minimum_loci=50, maf=0.05, max_missingness=0.3){
  
  out_matrix <- matrix(NA, 0, 12)
  colnames(out_matrix) <- c('genetic_group', 'site', 'Ar', 'Ho', 'He', 'uHe', 'Fis', 'uFis', 'loci','n', 'rAr','sd_rAr')
  
  for(group in unique(genetic_group_variable)){
    gt_group <- gt[which(genetic_group_variable==group),]
    
    not_missing_loci <- which(colMeans(is.na(gt_group))<=max_missingness)
    gt_group_missing <- gt_group[,not_missing_loci]
    
    loci_mafs <- get_minor_allele_frequencies(gt_group_missing)
    passing_maf_loci <- which(loci_mafs>=maf)
    gt_group_missing_maf <- gt_group_missing[,passing_maf_loci]
    
    if(ncol(gt_group_missing_maf)<minimum_loci){
      print(paste(group),'does not have enough loci (',ncol(gt_group_missing_maf),'loci: minimum is,',minimum_loci,')')
      print('Proceeding to next genetic group')
      next()
    }
    
    # Remove sites that have fewer samples than the minimum_n threshold
    site_freq <- table(site_variable)
    sites <- names(which(site_freq>=minimum_n))
    print(paste("Sites removed due to low n:",names(which(site_freq<minimum_n))))
    
    group_out_matrix <- matrix(NA, length(sites), 10)
    
    allele_counts <- list() 
    
    for(s in 1:length(sites)){
      site <- sites[s]
      gt_site <- gt_group_missing_maf[which(site_variable==site),]
      Ho <- calculate_Ho(gt_site) 
      Hes <- calculate_Hes(gt_site)
      He <- mean(Hes, na.rm=TRUE) 
      ns <- colSums(!is.na(gt_site))
      uHe <- calculate_uHe(ns, Hes) 
      Fis <- 1 - (Ho/He) 
      uFis <- 1 - (Ho/uHe) 
      loci <- ncol(gt_site)
      Ar <- calculate_Ar(gt_site)
      n <- nrow(gt_site)
      group_out_matrix[s,] <- c(group, site, 
                                round(Ar,3), round(Ho,3), round(He,3), round(uHe,3), round(Fis,3), round(uFis,3), 
                                loci, n)
    }
    
    ## calculate rarified allelic richness ##
    # rows are sites and columns are loci 
    count_matrix <- t(sapply(split(seq_along(site_variable), site_variable), 
                             function(idx) colSums(gt_group_missing_maf[idx, ], na.rm = TRUE))) %>% as.matrix()
    n_matrix <- t(sapply(split(seq_along(site_variable), site_variable), 
                         function(idx) colSums(!is.na(gt_group_missing_maf[idx, ]))))%>% as.matrix()
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
    
    # Calculate the mean allelic richness for each population, ignoring the first column and any NA values
    mn_ar<- rowMeans(all_ar, na.rm = TRUE) %>% round(.,3)
    
    # Calculate the standard deviation of allelic richness for each population, ignoring the first column and any NA values
    sd_ar <- apply(all_ar, 1, sd, na.rm = TRUE) %>% round(.,3)
    
    group_out_matrix <- cbind(group_out_matrix, mn_ar, sd_ar)
    
    out_matrix <- rbind(out_matrix, group_out_matrix)
    print(paste(group,"complete"))###
    
  }
  print('Process complete!')
  return(out_matrix%>% as.data.frame())
}