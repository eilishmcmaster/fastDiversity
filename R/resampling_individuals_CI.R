#' Calculate global bootstrapped He, Ho, and FIS by resampling individuals
#'
#' @description
#' Resampling individuals means randomly selecting individuals (rows) from the dataset with replacement. This keeps the same genetic markers (loci), but uses different combinations of individuals in each sample. It tests how your results might change depending on which individuals you sampled.
#' 
#' @param gt Genotype matrix where individuals are rows and loci are columns, coded as 
#'           0=aa, 1=aA, 2=AA.
#' @param boots Number of bootstrap replicates
#' @param resample_n Number of individuals to resample with replacement (defaults to total individuals)
#' @param CI_alpha Alpha level for confidence intervals (e.g., 0.05 for 95\% CI)
#' @return Named vector of confidence intervals for global FIS, HE, HO
#' @details
#' Bootstrap replicates are drawn with R's random number generator. Call
#' \code{set.seed()} before this function for reproducible confidence intervals.
#' @export
resampling_individuals_CI <- function(gt, boots = 100, resample_n = NULL, CI_alpha = 0.05) {
  
  if (is.null(resample_n)) resample_n <- nrow(gt)
  
  gt_mat <- as.matrix(gt)
  storage.mode(gt_mat) <- "double"
  n_ind <- nrow(gt_mat)
  
  # Bootstrap resampling of individuals
  resample_results <- lapply(seq_len(boots), function(i) {
    idx <- sample.int(n_ind, resample_n, replace = TRUE)
    resampled <- gt_mat[idx, , drop = FALSE]
    
    ho <- calculate_Ho_vec_cpp(resampled)
    he <- calculate_Hes_vec_cpp(resampled)
    
    fis <- 1 - (ho / he)
    fis[!is.finite(fis)] <- NA
    
    list(ho = ho, he = he, fis = fis)
  })
  
  # Combine results
  ho_mat  <- do.call(rbind, lapply(resample_results, `[[`, "ho"))
  he_mat  <- do.call(rbind, lapply(resample_results, `[[`, "he"))
  fis_mat <- do.call(rbind, lapply(resample_results, `[[`, "fis"))
  
  ci_probs <- c(CI_alpha / 2, 1 - CI_alpha / 2)
  
  global_ho_ci  <- c(quantile(rowMeans(ho_mat, na.rm = TRUE), probs = ci_probs), 
                     mean = mean(rowMeans(ho_mat, na.rm = TRUE), na.rm = TRUE))
  global_he_ci  <- c(quantile(rowMeans(he_mat, na.rm = TRUE), probs = ci_probs),
                     mean = mean(rowMeans(he_mat, na.rm = TRUE), na.rm = TRUE))
  global_fis_ci <- c(quantile(rowMeans(fis_mat, na.rm = TRUE), probs = ci_probs),
                     mean = mean(rowMeans(fis_mat, na.rm = TRUE), na.rm = TRUE))
  
  return(c('individuals_f' = round(global_fis_ci, 3),
           'individuals_he' = round(global_he_ci, 3),
           'individuals_ho' = round(global_ho_ci, 3)))
}
