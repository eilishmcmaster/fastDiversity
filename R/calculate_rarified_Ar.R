#' Calculate rarified allelic richness
#' 
#' This function calculates rarified allelic richness across genetic groups (sites) 
#' using rarefaction. Rarified allelic richness estimates the average number of 
#' alleles in a sample while accounting for differences in sample sizes.
#' 
#' References:
#' - Kalinowski, S. T. (2005). https://www.montana.edu/kalinowski/software/documents/2005_HP-Rare_MolecularEcologyNotes.pdf
#' - Keenan, K., et al. (2013). https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210x.12067
#' 
#' @param gt Genotype matrix where individuals are rows and loci are columns, coded as 
#'           0=aa, 1=aA, 2=AA.
#' @param site_variable Vector of group IDs (sites) for each individual (same length as 
#'                      number of rows in gt).
#' @return A matrix of rarefied allelic richness with one row per site and one column
#'         per locus. Each locus is rarefied to the smallest number of alleles sampled
#'         at that locus across sites. Loci with no genotypes at one or more sites are
#'         returned as NA.
#' @export
calculate_rarified_Ar <- function(gt, site_variable){
  ## calculate rarified allelic richness ##
  gt <- as.matrix(gt)
  site_rows <- split(seq_along(site_variable), site_variable)
  # rows are sites and columns are loci
  count_matrix <- t(vapply(site_rows, function(idx) colSums(gt[idx, , drop = FALSE], na.rm = TRUE),
                           numeric(ncol(gt)))) # alternate allele count
  n_matrix <- t(vapply(site_rows, function(idx) colSums(!is.na(gt[idx, , drop = FALSE])),
                       numeric(ncol(gt)))) # individuals genotyped
  if (ncol(gt) == 1) { # vapply drops to a vector for a single locus
    count_matrix <- t(count_matrix)
    n_matrix <- t(n_matrix)
  }
  a_matrix <- 2*n_matrix # alleles sampled per site and locus
  count_matrix2 <- a_matrix - count_matrix # reference allele count
  
  # rarefy each locus to the smallest number of alleles sampled at that locus
  g <- matrix(2 * apply(n_matrix, 2, min), nrow(n_matrix), ncol(n_matrix), byrow = TRUE)
  
  # probability that an allele is NOT seen in g draws without replacement,
  # on the log scale to avoid overflow in choose() for large sample sizes
  p_missing_alt <- exp(lchoose(a_matrix - count_matrix, g) - lchoose(a_matrix, g))
  p_missing_ref <- exp(lchoose(a_matrix - count_matrix2, g) - lchoose(a_matrix, g))
  all_ar <- (1 - p_missing_alt) + (1 - p_missing_ref) # expected number of alleles seen
  
  all_ar[, g[1, ] == 0] <- NA # locus not genotyped at every site
  dimnames(all_ar) <- list(names(site_rows), colnames(gt))
  return(all_ar)
}
