calculate_Ar <- function(gt){ #https://www.montana.edu/kalinowski/software/documents/2005_HP-Rare_MolecularEcologyNotes.pdf
  apc <- colSums(gt==0, na.rm=TRUE)+colSums(gt==2, na.rm=TRUE)+2*colSums(gt==1, na.rm=TRUE) # number of alleles at the locus
  ns <- colSums(!is.na(gt)) # number of individuals with the locus
  Ar <- mean(apc/ns, na.rm=TRUE) # average number of alleles per individual
  return(Ar)
}