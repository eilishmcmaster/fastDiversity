#' Identify Alleles Matching Genotype
#' 
#' This function identifies alleles at specified loci that match the genotype data (`gt`).
#' It checks the genotype matrix for each locus and retrieves corresponding alleles from 
#' a provided loci data frame (`loci`).
#' 
#' @param gt Genotype matrix where the first row represents loci and subsequent rows are genotype data
#'           coded as 0=aa, 1=aA, 2=AA.
#' @param loci Data frame containing locus information with columns 'loci', 'allele1', and 'allele2'.
#' @return Vector of alleles matching the genotype at each locus in `gt`.
#' @export
matcher <- function(gt, loci){
  df <- gt[-1]
  out <- vector()
  if(0 %in% df | 1 %in% df){ # allele 1 (a) is present in aa (0) and aA (1)
    out <- append(out, loci[loci==gt[1],2])
  }
  if(1 %in% df | 2 %in% df){ # allele 2 (A) is present in aA (1) and AA (2)
    out <- append(out,loci[loci==gt[1],3])
  }
  return(out)
}
