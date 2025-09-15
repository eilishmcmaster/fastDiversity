#' Test for Hardy-Weinberg Equilibrium (HWE) across multiple loci using HardyWeinberg
#'
#' This function performs Hardy-Weinberg equilibrium tests on multiple genetic loci
#' using the HWExact from the HardyWeinberg package. It returns the Fisher 
#' combined p-value of the LLR test, the mean U-score, and individual test results 
#' for each locus.
#'
#' @param gt A data frame or matrix where each column corresponds to genotypes at a single locus.
#'           Genotypes must be coded as 0 (AA), 1 (AB), or 2 (BB).
#'
#' @return A named vector containing:
#' \describe{
#'   \item{HW_fisher_p}{Fisher combined p-value from the HWExact test across loci.}
#'   \item{loci_HWX_results}{A data table with individual p-values and locus names.}
#' }
#'
#' @import data.table HardyWeinberg dplyr stats
#' @export
MC_HWEtest_of_loci <- function(gt) {
  
  # Helper function to generate genotype count matrix (AA, AB, BB) for one locus
  get_allele_freq_matrix <- function(geno_vector) {
    # Count each genotype (assumes values 0=AA, 1=AB, 2=BB)
    counts <- tabulate(factor(geno_vector, levels = c(0, 1, 2)), nbins = 3)
    # Create 2x2 symmetric matrix as expected by HWExact
    matrix(
      c(counts[1], counts[2],
        NA, counts[3]),
      nrow = 2, dimnames = list(c("A", "B"), c("A", "B"))
    )
  }
  
  # Apply Hardy-Weinberg test to each locus
  
  hwx_results <- lapply(seq_len(ncol(gt)), function(index) {
    geno_counts <- get_allele_freq_matrix(gt[[index]])
    if(isTRUE(sum(geno_counts, na.rm=TRUE)==0)){
      return(NULL)
    }
    result <- tryCatch({
      HardyWeinberg::HWExact(c(geno_counts[1,1], geno_counts[2,1], geno_counts[2,2]),
                             verbose=FALSE)
    }, error = function(e) {
      warning(paste("Error at locus", index, ":", e$message))
      return(NULL)
    })
    
    return(result)
  })
  
  # Assign locus names to the results
  names(hwx_results) <- names(gt)
  
  # Extract p-values from results
  HW_p <- sapply(hwx_results, function(x) if(!is.null(x$pval)) x$pval[1] else NA)
  
  # Fisher's method to combine p-values from the HW tests
  HW_X <- -2 * sum(log(HW_p), na.rm = TRUE)
  HW_df <- 2 * sum(!is.na(HW_p))
  fisher_combined_p <- pchisq(HW_X, HW_df, lower.tail = FALSE)
  
  # Return summary statistics and results
  return(list('HW_fisher_p' = round(fisher_combined_p, 6),
           'loci_HWX_results' = HW_p))
  
}
