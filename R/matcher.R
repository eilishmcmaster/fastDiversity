#' Identify Alleles Matching Genotype
#' 
#' This function identifies alleles at specified loci that match the genotype data (`gt`).
#' It checks the genotype matrix for each locus and retrieves corresponding alleles from 
#' a provided loci data frame (`loci`).
#' 
#' @param gt Genotype matrix where the first row represents loci and subsequent rows are genotype data (0, 1, 2).
#' @param loci Data frame containing locus information with columns 'loci', 'allele1', and 'allele2'.
#' @return Vector of alleles matching the genotype at each locus in `gt`.
#' @export
matcher <- function(df2, loci){
  df <- df2[-1]
  out <- vector()
  if(0 %in% df | 2 %in% df){
    out <- append(out, loci[loci==df2[1],2])
  }
  if(1 %in% df | 2 %in% df){
    out <- append(out,loci[loci==df2[1],3])
  }
  return(out)
}
