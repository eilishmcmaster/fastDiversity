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
#' @param site_variable (optional) Vector of site IDs for each individual. If supplied, expected heterozygosity will be calculated within sites, which will probably result in F being lower due to the Wahlund effect.
#' @param minimum_loci (optional) Minimum number of loci required to proceed with analysis for each 
#'                    genetic group (genetic groups with fewer loci will be skipped).
#' @param maf (optional) Minimum minor allele frequency filter.
#' @param max_missingness (optional) Maximum proportion of missing data allowed per locus (loci with 
#'                        higher missingness will be filtered out).
#' @return A named vector containing F for each individual in gt.
#' @examples
#' data(example_data)
#' # get inbreeding coefficients for individuals by genetic group
#' Fs <- individual_F(example_gt, example_meta$population)
#' @references 
#' Frankham, R., Ballou, J. D., & Briscoe, D. A. (2010). Introduction to Conservation Genetics (2nd ed.). Cambridge University Press. https://doi.org/10.1017/CBO9780511809002
#' @export
individual_F <- function(gt, 
                          genetic_group_variable = NULL, 
                          site_variable = NULL,
                          minimum_loci = 50, 
                          maf = 0.05, 
                          max_missingness = 0.3,
                          minimum_n = 3, 
                          table_out = FALSE) {
  
  if (!is.null(site_variable)) {
    site_variable <- as.character(site_variable)
  }
  if (!is.null(genetic_group_variable)) {
    genetic_group_variable <- as.character(genetic_group_variable)
  }
  
  # table to store outputs
  all_inds <- rownames(gt)
  results <- data.frame(individual = all_inds, 
                        F = NA_real_,
                        group = if (!is.null(genetic_group_variable)) genetic_group_variable else NA,
                        site = if (!is.null(site_variable)) site_variable else NA,
                        stringsAsFactors = FALSE)
  
  # subfunctions
  filter_loci_function <- function(gt){
    not_missing_loci <- which(colMeans(is.na(gt)) <= max_missingness)
    gt_missing <- gt[, not_missing_loci, drop=FALSE]
    
    loci_mafs <- get_minor_allele_frequencies(gt_missing)
    passing_maf_loci <- which(loci_mafs >= maf)
    gt_filtered <- gt_missing[, passing_maf_loci, drop=FALSE]
    return(gt_filtered)
  }
  
  
  calculate_f_function <- function(gt){
    if (is.null(gt) || !is.matrix(gt) && !is.data.frame(gt)) {
      message("Invalid input: gt must be a non-null matrix or data frame.")
      return(NA)
    }
    if (nrow(gt) < minimum_n) {
      message(paste0("Skipping ungrouped dataset - fewer than ", minimum_n, " samples (n = ", nrow(gt), ")"))
      return(rep(NA, nrow(gt)))
    }
    if (ncol(gt) < minimum_loci) {
      message(paste0("Dataset does not have enough loci (", ncol(gt), " loci: minimum is ", minimum_loci, ")"))
      return(rep(NA, nrow(gt)))
    }
      individual_heterozygosity <- rowSums((gt == 1), na.rm = TRUE) /rowSums(!is.na(gt))
      expected_heterozygosity <- mean(calculate_Hes(gt), na.rm = TRUE)
      
      F_values <- 1 - (individual_heterozygosity / expected_heterozygosity)
      return(F_values)
  }

  #### run for various scenarios ####
  
  # No population or site groups
  if(is.null(genetic_group_variable) & is.null(site_variable)){
    gt_filtered <- filter_loci_function(gt)
    F_values <- calculate_f_function(gt_filtered)
  }
  
  
  # no population but yes site groups
  if(is.null(genetic_group_variable) & !is.null(site_variable)){
    gt_filtered_0 <- filter_loci_function(gt)
    sites <- unique(site_variable)
    F_values <- c()
    for(site in sites){
      gt_filtered <- gt_filtered_0[which(site_variable==site),]
      F_values_site <- calculate_f_function(gt_filtered)
      F_values <- c(F_values, F_values_site)
    }
  }
  
  # yes population no site 
  if(!is.null(genetic_group_variable) & is.null(site_variable)){
    groups <- unique(genetic_group_variable)
    F_values <- c()
    for(group in groups){
      gt_group <- gt[which(genetic_group_variable==group),]
      gt_filtered <- filter_loci_function(gt_group)
      F_values_group <- calculate_f_function(gt_filtered)
      F_values <- c(F_values, F_values_group)
    }
  }
  
  # yes population yes site 
  if(!is.null(genetic_group_variable) & !is.null(site_variable)){
    groups <- unique(genetic_group_variable)
    F_values <- c()
    for(group in groups){
      gt_group <- gt[which(genetic_group_variable==group),]
      group_sites <- site_variable[which(genetic_group_variable==group)]
      gt_filtered <- filter_loci_function(gt_group)
      sites <- unique(group_sites)
      
      for(site in sites){
        gt_filtered_2 <- gt_filtered[which(group_sites==site),]
        F_values_site <- calculate_f_function(gt_filtered_2)
        F_values <- c(F_values, F_values_site)
      }
    }
  }
  
  # return(F_values)
  results$F <- F_values[results$individual]
  if(isFALSE(table_out)){
    rownames(results) <- NULL
    out_F <- results$F %>% as.vector()
    names(out_F) <- results$individual
    return(out_F)
  }
  if(isTRUE(table_out)){
    rownames(results) <- NULL
    return(results)
  }
}
