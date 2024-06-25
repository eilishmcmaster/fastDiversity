#' Calculate Unbiased Expected Heterozygosity (uHe)
#' 
#' This function calculates unbiased expected heterozygosity (uHe) based on the number 
#' of individuals (ns) and observed heterozygosity (Hes) for each genetic group or locus.
#' Unbiased expected heterozygosity corrects for sample size biases and provides an 
#' estimate of heterozygosity accounting for the number of alleles and individuals.
#' 
#' @param ns Vector containing the number of individuals sampled for each genetic group 
#'           or locus.
#' @param Hes Vector containing observed heterozygosity values for each genetic group 
#'            or locus.
#' @return Average unbiased expected heterozygosity (uHe) across genetic groups or loci.
#' @export
calculate_uHe <- function(ns, Hes){
  uHes <- c()
  for(i in 1:length(ns)){
    uHes[i] <- (2*ns[i]/(2*ns[i]-1))*Hes[i]
  }
  uHe <- mean(uHes,  na.rm=TRUE) 
  return(uHe)
}