#' Get Minor Allele Counts
#' 
#' This function calculates the minor allele counts for each locus based on genotype data.
#' The minor allele count is defined as the sum of the less common allele at each locus.
#' 
#' @param gt Genotype matrix where individuals are rows and loci are columns, coded as 
#'           0=aa, 1=aA, 2=AA.
#' @return Vector of minor allele counts for each locus.
#' @export
get_minor_allele_counts <- function(mat) {
  alt <- colSums(mat, na.rm = TRUE)
  n_called <- colSums(!is.na(mat))
  total_alleles <- 2 * n_called
  ref <- total_alleles - alt
  pmin(ref, alt)
}