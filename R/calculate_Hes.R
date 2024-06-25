calculate_Hes <- function(gt){
  alleles   <- 2*(colSums(!is.na(gt))) # total number of alleles for a locus (all samples in gt)
  alt_freq  <- colSums(gt,na.rm=TRUE) / (alleles) # get the allele frequencies (altcount data can be summed because 0=homo1, 1=het, 2=homo2)
  ref_freq  <- 1-alt_freq # get the alternative allele frequency
  
  Hes <- c()
  for(i in 1:ncol(gt)){
    Hes[i] <- 1-sum(ref_freq[i]^2+alt_freq[i]^2)
  }
  return(Hes)
}