#' Get Minor Allele Frequencies
#' 
#' This function calculates the minor allele frequencies for each locus based on genotype data.
#' The minor allele frequency is defined as the frequency of the less common allele at each locus.
#' 
#' @param gt Genotype matrix where individuals are rows and loci are columns, coded as 
#'           0=aa, 1=aA, 2=AA.
#' @return Vector of minor allele frequencies for each locus.
#' @export
get_minor_allele_frequencies <- function( gt ) {
  alleles   <- 2*(colSums(!is.na(gt))) # total number of alleles for a locus (all samples in gt)
  alt_freq  <- colSums(gt,na.rm=TRUE) / (alleles) # get the allele frequencies (altcount data can be summed because 0=homo1, 1=het, 2=homo2)
  ref_freq  <- 1-alt_freq # get the alternative allele frequency
  min_freq <- alt_freq
  for ( i in 1:ncol(gt) ) { # assign the minor allele to smaller value
    if ( isTRUE(alt_freq[i] > ref_freq[i] )) {
      min_freq[i] <- ref_freq[i]
    } 
  }
  return(min_freq) # return minor allele frequency for each locus in gt 
}