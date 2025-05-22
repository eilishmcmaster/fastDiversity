#' Calculate global bootstrapped He, Ho, and FIS
#' 
#' @param gt Genotype matrix where individuals are rows and loci are columns, coded as 
#'           0=aa, 1=aA, 2=AA.
#' @param fis_boots Number of bootstrap replicates
#' @param resample_n Number of individuals to resample from the site -- defaults to site n
#' @return Vector of observed heterozygosity values (Hes) for each locus.
#' @export
calculate_boot_stats <- function(gt, boots, resample_n){
  if(is.null(resample_n)){ # if no resample_n is supplied, use the number of samples in the site
    resample_n <- nrow(gt) 
  }
  if(is.null(fis_boots)){ # if no fis_boots is supplied, use 100
    fis_boots <- 100
  }
  
  resample_list <- lapply(1:fis_boots, function(i) sample(1:nrow(gt), resample_n, replace = TRUE))
  
  n_ind <- nrow(gt)
  loci <- names(gt)
  
  # initialize storage
  ho_list <- vector("list", fis_boots)
  he_list <- vector("list", fis_boots)
  fis_list <- vector("list", fis_boots)
  
  loci_Ho <- function(x) {
    if (all(is.na(x))) {
      return(NA)
    } else {
      return(mean(x == 1, na.rm = TRUE))
    }
  }
  
  # bootstrapping
  for (i in seq_len(fis_boots)) {
    resample_idx <- resample_list[[i]] #sample(n_ind, n_ind, replace = TRUE)
    resampled_dt <- gt[resample_idx]
    
    ho <- resampled_dt[, lapply(.SD, loci_Ho)] # heterozygosity per locus
    he <- calculate_Hes(resampled_dt) # expected heterozygosity of a column
    
    fis <- 1 - (as.numeric(ho) / as.numeric(he))
    fis[is.nan(fis)] <- NA
    
    ho_list[[i]] <- unlist(ho)
    he_list[[i]] <- unlist(he)
    fis_list[[i]] <- fis
  }
  
  # convert lists to matrices
  ho_mat <- do.call(rbind, ho_list)
  he_mat <- do.call(rbind, he_list)
  fis_mat <- do.call(rbind, fis_list)
  
  # convert to data.tables for CI calculations
  ho_dt <- as.data.table(ho_mat) %>% setnames(., loci)
  he_dt <- as.data.table(he_mat) %>% setnames(., loci)
  fis_dt <- as.data.table(fis_mat) %>% setnames(., loci)
  
  ### locus fis 
  # get fis per locus x boots
  # calculate ci for each locus across the boots (results in ci x loci)
  ho_ci <- ho_dt[, lapply(.SD, function(x) quantile(x, probs = c(fis_alpha/2, 1 - fis_alpha/2), na.rm = TRUE))] # CI of each locus
  he_ci <- he_dt[, lapply(.SD, function(x) quantile(x, probs = c(fis_alpha/2, 1 - fis_alpha/2), na.rm = TRUE))]
  fis_ci <- fis_dt[, lapply(.SD, function(x) quantile(x, probs = c(fis_alpha/2, 1 - fis_alpha/2), na.rm = TRUE))]
  
  ### global fis
  # get fis per locus x boots
  # get the mean fis of each locus
  # get cis across all bootstraps
  global_fis_ci <- quantile(rowMeans(fis_dt, na.rm = TRUE), c(fis_alpha/2, 1 - fis_alpha/2)) # CI of mean Fis across boots
  global_he_ci <- quantile(rowMeans(he_dt, na.rm = TRUE), c(fis_alpha/2, 1 - fis_alpha/2)) # CI of mean He across boots
  global_ho_ci <- quantile(rowMeans(ho_dt, na.rm = TRUE), c(fis_alpha/2, 1 - fis_alpha/2)) # CI of mean Ho across boots
  
  return(c('global_fis_ci' = round(global_fis_ci, 3), 
           'global_he_ci' = round(global_he_ci, 3), 
           'global_ho_ci' = round(global_ho_ci, 3)))
  
}
