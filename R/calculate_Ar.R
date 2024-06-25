#' Calculate Allelic Richness
#' 
#' This function calculates allelic richness (Ar) for each locus based on genotype data.
#' Allelic richness is defined as the average number of alleles per individual across all loci.
#' 
#' Reference:
#' - Kalinowski, S. T. (2005). https://www.montana.edu/kalinowski/software/documents/2005_HP-Rare_MolecularEcologyNotes.pdf
#' 
#' @param gt Genotype matrix where individuals are rows and loci are columns, coded as 
#'           0=aa, 1=aA, 2=AA.
#' @return Average allelic richness across loci.
#' @export
calculate_Ar <- function(gt){ 
  apc <- colSums(gt==0, na.rm=TRUE)+colSums(gt==2, na.rm=TRUE)+2*colSums(gt==1, na.rm=TRUE) # number of alleles at the locus
  ns <- colSums(!is.na(gt)) # number of individuals with the locus
  Ar <- mean(apc/ns, na.rm=TRUE) # average number of alleles per individual
  return(Ar)
}