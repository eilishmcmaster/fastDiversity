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