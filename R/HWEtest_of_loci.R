#' Test for Hardy-Weinberg Equilibrium (HWE) across multiple loci using HWxtest
#'
#' This function performs Hardy-Weinberg equilibrium tests on multiple genetic loci
#' using the hwx.test from the HWxtest package. It returns the Fisher 
#' combined p-value of the LLR test, the mean U-score, and individual test results 
#' for each locus.
#'
#' @param gt A data frame or matrix where each column corresponds to genotypes at a single locus.
#'           Genotypes must be coded as 0 (AA), 1 (AB), or 2 (BB).
#'
#' @return A named vector containing:
#' \describe{
#'   \item{LLR_fisher_p}{Fisher combined p-value from the log-likelihood ratio (LLR) test across loci.}
#'   \item{mean_U_score}{Mean of the U-scores across all loci.}
#'   \item{loci_HWX_results}{A data table with individual LLR p-values, U-scores, and locus names.}
#' }
#'
#' @importFrom data.table data.table
#' @importFrom HWxtest hwx.test
#' @examples
#' @export
HWEtest_of_loci <- function(gt) {
  
  # Helper function to generate genotype count matrix (AA, AB, BB) for one locus
  get_allele_freq_matrix <- function(geno_vector) {
    # Count each genotype (assumes values 0=AA, 1=AB, 2=BB)
    counts <- tabulate(factor(geno_vector, levels = c(0, 1, 2)), nbins = 3)
    # Create 2x2 symmetric matrix as expected by hwx.test
    matrix(
      c(counts[1], counts[2],
        NA, counts[3]),
      nrow = 2, dimnames = list(c("A", "B"), c("A", "B"))
    )
  }
  
  # Apply Hardy-Weinberg test to each locus
  hwx_results <- lapply(seq_len(ncol(gt)), function(index) {
    geno_counts <- get_allele_freq_matrix(gt[[index]])
    result <- tryCatch({
      HWxtest::hwx.test(geno_counts)
    }, error = function(e) {
      warning(paste("Error at locus", index, ":", e$message))
      return(NULL)
    })
    return(result)
  })
  
  # Assign locus names to the results
  names(hwx_results) <- names(gt)
  
  # Extract LLR p-values from results
  LLR_p <- lapply(hwx_results, function(x) x[["Pvalues"]][["LLR"]]) %>% unlist
  
  # Extract U-scores (heterozygote or homozygote excess)
  U_score <- lapply(hwx_results, function(x) x[["observed"]][["U"]]) %>% unlist
  
  # Fisher's method to combine p-values from the LLR tests
  LLR_X <- -2 * sum(log(LLR_p))
  LLR_df <- 2 * length(LLR_p)
  fisher_combined_p <- pchisq(LLR_X, LLR_df, lower.tail = FALSE)
  
  # Create summary data table for each locus
  loci_results_dt <- data.table(LLR_p, U_score, locus = names(gt))
  
  # Return summary statistics and results
  return(c('LLR_fisher_p' = fisher_combined_p,
           'mean_U_score' = mean(U_score),
           'loci_HWX_results' = list(loci_results_dt)))
}
