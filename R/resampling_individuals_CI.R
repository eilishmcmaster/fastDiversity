#' Calculate global bootstrapped He, Ho, and FIS (parallelized) by resampling individuals (as in diveRsity::basicstats)
#' 
#' @description
#' Resampling individuals means randomly selecting individuals (rows) from the dataset with replacement. This keeps the same genetic markers (loci), but uses different combinations of individuals in each sample. It tests how your results might change depending on which individuals you sampled.
#' Resampling loci means randomly selecting genetic markers (columns) with replacement. This keeps the same individuals, but uses different combinations of loci in each sample. It tests how your results might change depending on which loci were included.
#' 
#' @param gt Genotype matrix where individuals are rows and loci are columns, coded as 
#'           0=aa, 1=aA, 2=AA.
#' @param boots Number of bootstrap replicates
#' @param resample_n Number of individuals to resample (with replacement) -- defaults to site n
#' @param CI_alpha Alpha level for confidence intervals (e.g., 0.05 for 95% CI)
#' @return Named vector of confidence intervals for global FIS, HE, HO
#' @import data.table dplyr future.apply
#' @export
resampling_individuals_CI <- function(gt, boots, resample_n, CI_alpha = 0.05) {
  
  # If no resample_n is supplied, use the number of samples in the site
  if (is.null(resample_n)) resample_n <- nrow(gt)
  if (is.null(boots)) boots <- 100
  
  # Function to compute observed heterozygosity at each locus:
  # Ho = proportion of heterozygotes (coded as 1)
  loci_Ho <- function(x) {
    if (all(is.na(x))) {
      return(NA)
    } else {
      return(mean(x == 1, na.rm = TRUE))
    }
  }
  
  # Convert genotype matrix to data.table for fast row/column manipulation
  gt_dt <- as.data.table(gt)
  loci <- names(gt)
  
  # Setup parallel processing using 'future' framework
  # multisession is safe across OS platforms (use multicore for Linux if desired)
  future::plan(future::multisession)
  
  # Bootstrap loop:
  # For each replicate:
  # 1. Resample individuals with replacement
  # 2. Calculate Ho (observed heterozygosity)
  # 3. Calculate He (expected heterozygosity)
  # 4. Calculate FIS = 1 - (Ho / He)
  # Results are stored in a list of length boots, with each element a list of Ho, He, FIS
  # system.time({
  resample_results <- future_lapply(1:boots, function(i) {
    resample_idx <- sample(nrow(gt_dt), resample_n, replace = TRUE)
    resampled_dt <- gt_dt[resample_idx]
    
    ho <- resampled_dt[, vapply(.SD, loci_Ho, numeric(1))] # Ho per locus
    he <- calculate_Hes(resampled_dt)                      # He per locus (custom function)
    
    fis <- 1 - (ho / he)                                   # FIS per locus
    fis[is.nan(fis)] <- NA                                 # Set NaNs to NA
    
    list(ho = ho, he = he, fis = fis)
  }, future.seed = 14)  # Ensures reproducibility of bootstrapping
  # })
  
  # Extract Ho, He, and FIS values from each bootstrap iteration and bind into matrices
  # Each matrix will be (boots x loci)
  ho_mat  <- do.call(rbind, lapply(resample_results, `[[`, "ho"))
  he_mat  <- do.call(rbind, lapply(resample_results, `[[`, "he"))
  fis_mat <- do.call(rbind, lapply(resample_results, `[[`, "fis"))
  
  # Define lower and upper quantile bounds for confidence intervals
  ci_probs <- c(CI_alpha / 2, 1 - CI_alpha / 2)
  
  # ### Locus-level confidence intervals:
  # # For each locus, calculate the bootstrap confidence interval for Ho, He, and FIS
  # ho_ci  <- apply(ho_mat,  2, quantile, probs = ci_probs, na.rm = TRUE)
  # he_ci  <- apply(he_mat,  2, quantile, probs = ci_probs, na.rm = TRUE)
  # fis_ci <- apply(fis_mat, 2, quantile, probs = ci_probs, na.rm = TRUE)
  
  ### Global (mean across loci) confidence intervals:
  # 1. Take the mean across loci for each bootstrap replicate
  # 2. Calculate the CI across all replicates for each metric
  global_ho_ci  <- c(quantile(rowMeans(ho_mat,  na.rm = TRUE), probs = ci_probs), 'mean'=mean(rowMeans(ho_mat,  na.rm = TRUE)))
  global_he_ci  <- c(quantile(rowMeans(he_mat,  na.rm = TRUE), probs = ci_probs), 'mean'=mean(rowMeans(he_mat,  na.rm = TRUE)))
  global_fis_ci <- c(quantile(rowMeans(fis_mat, na.rm = TRUE), probs = ci_probs), 'mean'=mean(rowMeans(fis_mat,  na.rm = TRUE)))
  
  # Return named vector of rounded global confidence intervals
  return(c('individuals_f' = round(global_fis_ci, 3), 
           'individuals_he'  = round(global_he_ci, 3), 
           'individuals_ho'  = round(global_ho_ci, 3)))
}
