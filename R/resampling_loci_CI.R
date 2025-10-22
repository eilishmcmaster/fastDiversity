#' Calculate global bootstrapped He, Ho, and FIS (parallelized) by resampling loci (as in hierfstat::boot.ppfis or arlequin)
#'
#' @description
#' Resampling loci means randomly selecting genetic markers (columns) with replacement. This keeps the same individuals, but uses different combinations of loci in each sample. It tests how your results might change depending on which loci were included.
#' 
#' @param gt Genotype matrix where individuals are rows and loci are columns, coded as 
#'           0=aa, 1=aA, 2=AA.
#' @param boots Number of bootstrap replicates
#' @param resample_n Number of loci to resample with replacement -- defaults to total loci
#' @param CI_alpha Alpha level for confidence intervals (e.g., 0.05 for 95% CI)
#' @return Named vector of confidence intervals for global FIS, HE, HO
#' @import data.table dplyr future.apply Rcpp
#' @export
resampling_loci_CI <- function(gt, boots = 100, resample_n = NULL, CI_alpha = 0.05) {
  
  if (is.null(resample_n)) resample_n <- ncol(gt)
  
  ###
  
  cppFunction('
    NumericVector calculate_Ho_vec_cpp(NumericMatrix gt_mat) {
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
  
  cppFunction('
    NumericVector calculate_Hes_vec_cpp(NumericMatrix gt_mat) {
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
          double alt_freq = (double)alt_count / allele_count;
          he[j] = 2 * alt_freq * (1 - alt_freq);
        } else {
          he[j] = NA_REAL;
        }
      }
      return he;
    }
    ')
  
  ###
  gt_mat <- as.matrix(gt)
  storage.mode(gt_mat) <- "integer"  # force integer type
  loci_count <- ncol(gt_mat)
  
  # Bootstrap loop:
  # For each replicate:
  # 1. Resample loci with replacement
  # 2. Calculate Ho (observed heterozygosity)
  # 3. Calculate He (expected heterozygosity) 
  # 4. Calculate FIS = 1 - (Ho / He)
  # Results are stored in a list of length boots, with each element a list of Ho, He, FIS
  set.seed(14)
  resample_results <- lapply(seq_len(boots), function(i) {
    resample_idx <- sample.int(loci_count, resample_n, replace = TRUE)
    resampled <- gt_mat[, resample_idx, drop = FALSE]
    
    ho  <- calculate_Ho_vec_cpp(resampled)
    he  <- calculate_Hes_vec_cpp(resampled)
    
    fis <- 1 - (ho / he)
    fis[!is.finite(fis)] <- NA
    
    list(ho = ho, he = he, fis = fis)
  })
  
  # Combine results
  ho_mat  <- do.call(rbind, lapply(resample_results, `[[`, "ho"))
  he_mat  <- do.call(rbind, lapply(resample_results, `[[`, "he"))
  fis_mat <- do.call(rbind, lapply(resample_results, `[[`, "fis"))
  
  ci_probs <- c(CI_alpha / 2, 1 - CI_alpha / 2)
  
  global_ho_ci  <- c(quantile(rowMeans(ho_mat,  na.rm = TRUE), probs = ci_probs), 'mean'=mean(rowMeans(ho_mat,  na.rm = TRUE)))
  global_he_ci  <- c(quantile(rowMeans(he_mat,  na.rm = TRUE), probs = ci_probs), 'mean'=mean(rowMeans(he_mat,  na.rm = TRUE)))
  global_fis_ci <- c(quantile(rowMeans(fis_mat, na.rm = TRUE), probs = ci_probs), 'mean'=mean(rowMeans(fis_mat, na.rm = TRUE)))
  
  return(c('loci_f' = round(global_fis_ci, 3), 
           'loci_he'  = round(global_he_ci, 3), 
           'loci_ho'  = round(global_ho_ci, 3)))
}