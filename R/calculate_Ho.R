#' Calculate Observed Heterozygosity (Ho)
#' 
#' This function calculates observed heterozygosity (Ho) based on genotype data.
#' Heterozygosity (Ho) is defined as the proportion of heterozygous individuals 
#' among all individuals with genotype data.
#' 
#' @param gt Genotype matrix where individuals are rows and loci are columns, coded as 
#'           0=aa, 1=aA, 2=AA.
#' @return Observed heterozygosity (Ho) across all loci.
#' @export
calculate_Ho <- function(gt){
  heteros <- sum(gt==1, na.rm=TRUE)
  total <- sum(!is.na(gt))
  Ho <- heteros/total
  return(Ho)
}