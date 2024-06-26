#' Calculate individual inbreeding coefficient (Wright's F)
#' 
#' This function calculates Wright's F (individual inbreeding coefficient), the probability 
#' that two alleles at a locus in an individual are identical by descent. Calculation from
#' Frankham, Ballou, & Briscoe (2010, p 266) F= 1-(Ho/He).
#' 
#' @param gt Genotype matrix where individuals are rows and loci are columns, coded as 
#'           0=aa, 1=aA, 2=AA.
#' @param genetic_group_variable (optional) Vector of genetic group IDs for each individual (same 
#'                               length as number of rows in gt). Can represent species, 
#'                               subspecies, or other biologically relevant structure groups. If not
#'                               included, the whole gt will be treated as one group.
#' @param minimum_loci (optional) Minimum number of loci required to proceed with analysis for each 
#'                    genetic group (genetic groups with fewer loci will be skipped).
#' @param maf (optional) Minimum minor allele frequency filter.
#' @param max_missingness (optional) Maximum proportion of missing data allowed per locus (loci with 
#'                        higher missingness will be filtered out).
#' @return A named vector containing F for each individual in gt.
#' @examples
#' # get inbreeding coefficients for individuals by genetic group
#' Fs <- individual_F(gt, genetic_group_variable)
#' @references 
#' Frankham, R., Ballou, J. D., & Briscoe, D. A. (2010). Introduction to Conservation Genetics (2nd ed.). Cambridge University Press. https://doi.org/10.1017/CBO9780511809002
#' @export
individual_F <- function(gt, genetic_group_variable=NULL, minimum_loci=50, maf=0.05, max_missingness=0.3){
  if(is.null(genetic_group_variable)){
    not_missing_loci <- which(colMeans(is.na(gt_group))<=max_missingness)
    gt_group_missing <- gt_group[,not_missing_loci]
    
    loci_mafs <- get_minor_allele_frequencies(gt_group_missing)
    passing_maf_loci <- which(loci_mafs>=maf)
    gt_group_missing_maf <- gt_group_missing[,passing_maf_loci]
    
    if(ncol(gt_group_missing_maf)<minimum_loci){
      print(paste(group),'does not have enough loci (',ncol(gt_group_missing_maf),'loci: minimum is,',minimum_loci,')')
      print('Proceeding to next genetic group')
      next()
    }
    
    individual_heterozygosity <- rowSums((gt_group_missing_maf==1), na.rm=TRUE) /
      rowSums(!is.na(gt_group_missing_maf)) # number of hets per individual
    expected_heterozygosity <- mean(calculate_Hes(gt_group_missing_maf))
    
    F <- 1- (individual_heterozygosity/expected_heterozygosity)
  }
  else{
    F <- c()
    for(group in unique(genetic_group_variable)){
      gt_group <- gt[which(genetic_group_variable==group),]
      
      not_missing_loci <- which(colMeans(is.na(gt_group))<=max_missingness)
      gt_group_missing <- gt_group[,not_missing_loci]
      
      loci_mafs <- get_minor_allele_frequencies(gt_group_missing)
      passing_maf_loci <- which(loci_mafs>=maf)
      gt_group_missing_maf <- gt_group_missing[,passing_maf_loci]
      
      if(ncol(gt_group_missing_maf)<minimum_loci){
        print(paste(group),'does not have enough loci (',ncol(gt_group_missing_maf),'loci: minimum is,',minimum_loci,')')
        print('Proceeding to next genetic group')
        next()
      }
      individual_heterozygosity <- rowSums((gt_group_missing_maf==1), na.rm=TRUE) / rowSums(!is.na(gt_group_missing_maf)) # number of hets per individual
      expected_heterozygosity <- mean(calculate_Hes(gt_group_missing_maf))
      
      F_group <- 1- (individual_heterozygosity/expected_heterozygosity)
      F <- c(F, F_group)
    }
  }
  return(F)
}