#' Calculate global bootstrapped He, Ho, and FIS (parallelized) by resampling individuals
#'
#' @description
#' Resampling individuals means randomly selecting individuals (rows) from the dataset with replacement. This keeps the same genetic markers (loci), but uses different combinations of individuals in each sample. It tests how your results might change depending on which individuals you sampled.
#' 
#' @param gt Genotype matrix where individuals are rows and loci are columns, coded as 
#'           0=aa, 1=aA, 2=AA.
#' @param boots Number of bootstrap replicates
#' @param resample_n Number of individuals to resample with replacement (defaults to total individuals)
#' @param CI_alpha Alpha level for confidence intervals (e.g., 0.05 for 95% CI)
#' @return Named vector of confidence intervals for global FIS, HE, HO
#' @import data.table dplyr future.apply Rcpp
#' @export
resampling_individuals_CI <- function(gt, boots = 100, resample_n = NULL, CI_alpha = 0.05) {
  
  if (is.null(resample_n)) resample_n <- nrow(gt)
  
  # Rcpp functions
  Rcpp::cppFunction('
    NumericVector calculate_Ho_cpp(NumericMatrix gt_mat) {
      int n_col = gt_mat.ncol();
      int n_row = gt_mat.nrow();
      NumericVector ho(n_col);
      
      for (int j = 0; j < n_col; ++j) {
        int heteros = 0;
        int total = 0;
        for (int i = 0; i < n_row; ++i) {
          if (!NumericVector::is_na(gt_mat(i,j))) {
            total++;
            if (gt_mat(i,j) == 1) heteros++;
          }
        }
        ho[j] = total > 0 ? (double)heteros / total : NA_REAL;
      }
      return ho;
    }
  ')
  
  Rcpp::cppFunction('
    NumericVector calculate_He_cpp(NumericMatrix gt_mat) {
      int n_col = gt_mat.ncol();
      int n_row = gt_mat.nrow();
      NumericVector he(n_col);
      
      for (int j = 0; j < n_col; ++j) {
        int allele_count = 0;
        int alt_count = 0;
        for (int i = 0; i < n_row; ++i) {
          if (!NumericVector::is_na(gt_mat(i,j))) {
            allele_count += 2;
            alt_count += gt_mat(i,j);
          }
        }
        if (allele_count > 0) {
          double p = (double)alt_count / allele_count;
          he[j] = 2 * p * (1 - p);
        } else {
          he[j] = NA_REAL;
        }
      }
      return he;
    }
  ')
  
  # Ensure matrix format
  gt_mat <- as.matrix(gt)
  storage.mode(gt_mat) <- "integer"
  n_ind <- nrow(gt_mat)
  
  # Bootstrap resampling of individuals
  set.seed(14)
  resample_results <- lapply(seq_len(boots), function(i) {
    idx <- sample.int(n_ind, resample_n, replace = TRUE)
    resampled <- gt_mat[idx, , drop = FALSE]
    
    ho <- calculate_Ho_cpp(resampled)
    he <- calculate_He_cpp(resampled)
    
    fis <- 1 - (ho / he)
    fis[!is.finite(fis)] <- NA
    
    list(ho = ho, he = he, fis = fis)
  })
  
  # Combine results
  ho_mat  <- do.call(rbind, lapply(resample_results, `[[`, "ho"))
  he_mat  <- do.call(rbind, lapply(resample_results, `[[`, "he"))
  fis_mat <- do.call(rbind, lapply(resample_results, `[[`, "fis"))
  
  ci_probs <- c(CI_alpha / 2, 1 - CI_alpha / 2)
  
  global_ho_ci  <- c(quantile(rowMeans(ho_mat, na.rm = TRUE), probs = ci_probs), mean = mean(rowMeans(ho_mat, na.rm = TRUE)))
  global_he_ci  <- c(quantile(rowMeans(he_mat, na.rm = TRUE), probs = ci_probs), mean = mean(rowMeans(he_mat, na.rm = TRUE)))
  global_fis_ci <- c(quantile(rowMeans(fis_mat, na.rm = TRUE), probs = ci_probs), mean = mean(rowMeans(fis_mat, na.rm = TRUE)))
  
  return(c('individuals_f' = round(global_fis_ci, 3),
           'individuals_he' = round(global_he_ci, 3),
           'individuals_ho' = round(global_ho_ci, 3)))
}
