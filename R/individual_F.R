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
                         minimum_n = 3) {
  
  all_inds <- rownames(gt)
  results <- data.frame(individual = all_inds, 
                        F = NA_real_,
                        group = if (!is.null(genetic_group_variable)) genetic_group_variable else NA,
                        site = if (!is.null(site_variable)) site_variable else NA,
                        stringsAsFactors = FALSE)
  
  if (is.null(genetic_group_variable)) {
    # No grouping
    if (nrow(gt) < minimum_n) {
      message(paste("Skipping ungrouped dataset - fewer than", minimum_n, "samples (n =", nrow(gt), ")"))
      return(results)
    }
    
    not_missing_loci <- which(colMeans(is.na(gt)) <= max_missingness)
    gt_missing <- gt[, not_missing_loci, drop=FALSE]
    
    loci_mafs <- get_minor_allele_frequencies(gt_missing)
    passing_maf_loci <- which(loci_mafs >= maf)
    gt_filtered <- gt_missing[, passing_maf_loci, drop=FALSE]
    
    if (ncol(gt_filtered) < minimum_loci) {
      message(paste("Dataset does not have enough loci (", ncol(gt_filtered), 
                    " loci: minimum is ", minimum_loci, ")"))
      return(results)
    }
    
    individual_heterozygosity <- rowSums((gt_filtered == 1), na.rm = TRUE) /
      rowSums(!is.na(gt_filtered))
    expected_heterozygosity <- mean(calculate_Hes(gt_filtered))
    
    results$F <- 1 - (individual_heterozygosity / expected_heterozygosity)
    return(results)
  }
  
  # Grouped analysis
  for (group in unique(genetic_group_variable)) {
    group_idx <- which(genetic_group_variable == group)
    
    if (length(group_idx) < minimum_n) {
      message(paste("Skipping group", group, "- fewer than", minimum_n, "samples (n =", length(group_idx), ")"))
      next()
    }
    
    gt_group <- gt[group_idx, , drop=FALSE]
    
    if (!is.null(site_variable)) {
      # Group + Site nested loop
      sites_in_group <- unique(site_variable[group_idx])
      
      for (site in sites_in_group) {
        site_idx <- group_idx[which(site_variable[group_idx] == site)]
        
        if (length(site_idx) < minimum_n) {
          message(paste("Skipping group", group, "site", site, 
                        "- fewer than", minimum_n, "samples (n =", length(site_idx), ")"))
          next()
        }
        
        gt_site <- gt[site_idx, , drop=FALSE]
        
        not_missing_loci <- which(colMeans(is.na(gt_site)) <= max_missingness)
        gt_missing <- gt_site[, not_missing_loci, drop=FALSE]
        
        loci_mafs <- get_minor_allele_frequencies(gt_missing)
        passing_maf_loci <- which(loci_mafs >= maf)
        gt_filtered <- gt_missing[, passing_maf_loci, drop=FALSE]
        
        if (ncol(gt_filtered) < minimum_loci) {
          message(paste("Group", group, "Site", site, 
                        "does not have enough loci (", ncol(gt_filtered), 
                        " loci: minimum is ", minimum_loci, ")"))
          next()
        }
        
        individual_heterozygosity <- rowSums((gt_filtered == 1), na.rm = TRUE) /
          rowSums(!is.na(gt_filtered))
        expected_heterozygosity <- mean(calculate_Hes(gt_filtered))
        
        F_values <- 1 - (individual_heterozygosity / expected_heterozygosity)
        results$F[site_idx] <- F_values
      }
      
    } else {
      # Only group level
      not_missing_loci <- which(colMeans(is.na(gt_group)) <= max_missingness)
      gt_missing <- gt_group[, not_missing_loci, drop=FALSE]
      
      loci_mafs <- get_minor_allele_frequencies(gt_missing)
      passing_maf_loci <- which(loci_mafs >= maf)
      gt_filtered <- gt_missing[, passing_maf_loci, drop=FALSE]
      
      if (ncol(gt_filtered) < minimum_loci) {
        message(paste("Group", group, "does not have enough loci (", 
                      ncol(gt_filtered), " loci: minimum is ", minimum_loci, ")"))
        next()
      }
      
      individual_heterozygosity <- rowSums((gt_filtered == 1), na.rm = TRUE) /
        rowSums(!is.na(gt_filtered))
      expected_heterozygosity <- mean(calculate_Hes(gt_filtered))
      
      F_values <- 1 - (individual_heterozygosity / expected_heterozygosity)
      results$F[group_idx] <- F_values
    }
  }
  
  rownames(results) <- NULL
  return(results$F)
}