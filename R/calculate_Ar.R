#' Calculate Allelic Richness (unrarefied)
#'
#' This function calculates allelic richness (Ar) as the number of distinct alleles observed
#' at each locus (1 or 2 for biallelic SNPs), averaged across loci. This value is not
#' corrected for sample size, so it should only be compared between samples of similar size.
#' Use \code{\link{calculate_rarified_Ar}} for allelic richness standardised to a common
#' sample size, as reported by \code{\link{faststats}}.
#'
#' @param gt Genotype matrix where individuals are rows and loci are columns, coded as
#'           0=aa, 1=aA, 2=AA.
#' @return Mean number of alleles observed per locus. Loci with no genotypes are ignored.
#' @seealso \code{\link{calculate_rarified_Ar}}
#' @export
calculate_Ar <- function(gt){
  allele1_present <- colSums(gt == 0 | gt == 1, na.rm = TRUE) > 0 # a seen in aa or aA
  allele2_present <- colSums(gt == 1 | gt == 2, na.rm = TRUE) > 0 # A seen in aA or AA
  genotyped <- colSums(!is.na(gt)) > 0
  Ar <- mean((allele1_present + allele2_present)[genotyped])
  return(Ar)
}
