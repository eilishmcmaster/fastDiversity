#' Calculate pi, theta, and D statistics
#' 
#' @description
#' This function calculates pi, Wu and Watterson's theta, and Tajima's D
#' for specified sites within genetic groups following the methods in Pixy,
#' which accounts for missing data to produce unbiased estimates. For further
#' details see https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13326 and
#' https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.14104 . Data is
#' filtered based on genetic group (minimum number of loci, minimum minor 
#' allele frequency, and maximum missingness per locus) and minimum number of
#' individuals per site.
#' 
#' @param gt Genotype matrix where individuals are rows and loci are columns, coded as 
#'   0=aa, 1=aA, 2=AA.
#' @param site_variable Vector of group IDs (sites) for each individual (same length as 
#'   number of rows in gt).
#' @param genetic_group_variable Vector of genetic group IDs for each individual (same 
#'   length as number of rows in gt). Can represent species, subspecies, or other 
#'   biologically relevant structure groups.
#' @param minimum_n Minimum number of individuals required per site (sites with fewer 
#'   individuals will be filtered out).
#'  @param n Number of individuals to use for calculations per site. If left as NULL, all individuals are used.
#' @param minimum_loci Minimum number of loci required to proceed with analysis for each 
#'   genetic group (genetic groups with fewer loci will be skipped).
#' @param maf Minimum minor allele frequency filter.
#' @param allele_count_min Minimum minor allele count to keep a locus, can be used to exclude singletons. 
#' @param max_missingness Maximum proportion of missing data allowed per locus (loci with 
#'   higher missingness will be filtered out).
#' 
#' @return A data frame with the following columns:
#'   \item{group}{Genetic group ID (e.g., species, subspecies).}
#'   \item{site}{Site ID corresponding to the sample location.}
#'   \item{loci}{Number of loci used in the analysis for that group-site combination.}
#'   \item{n}{Number of individuals sampled at the site for the given group.}
#' @import data.table HardyWeinberg dplyr stats
#' @seealso [resampling_individuals_CI(), resampling_loci_CI()]
#' 
#' 
#' @export
pi_theta_d_stats <- function (gt, genetic_group_variable, site_variable, 
                       n = NULL, minimum_n = 5,  
                       minimum_loci = 5, max_missingness = 0.1, 
                       # presample_genetic_groups = NULL,
                       maf = 0, allele_count_min = 0) {
  
  gt <- as.data.table(gt, keep.rownames = FALSE)
  setnames(gt, make.names(colnames(gt), unique=TRUE))
  
  out_list <- list() # make empty list to put data in 
  
  genetic_group_freq <- table(genetic_group_variable) # get number of samples per genetic group
  genetic_groups <- names(which(genetic_group_freq >= minimum_n)) # get the genetic groups with n>=minimum
  
  for (group in genetic_groups) { # for each genetic group
    gt_group <- gt[which(genetic_group_variable == group), ] # filter samples
    site_variable_group <- site_variable[which(genetic_group_variable ==  group)] # get vector of sites for samples in the group
    #
    site_freq <- table(site_variable_group) # get number of samples per site
    sites <- names(which(site_freq >= minimum_n)) # remove sites less than threshold n
    
    if (length(sites) == 0) {
      message(sprintf("Skipping group '%s': no sites meet n ≥ %d.", 
                      group, minimum_n)) 
      next  }
    
    # gt_group_missing <- gt_group #TRIAL
    
    not_missing_loci <- which(colMeans(is.na(gt_group)) <=  max_missingness)  # determine which loci pass missingness threshold
    gt_group_missing <- gt_group[, ..not_missing_loci] # filter to remove high missing loci
    if (ncol(gt_group_missing) == 0) next
    
    if(maf>0){
      loci_mafs <- get_minor_allele_frequencies(gt_group_missing) # get minor allele frequencies of each locus
      passing_maf_loci <- which(loci_mafs >= maf) # get indices
      gt_group_missing <- gt_group_missing[, ..passing_maf_loci] # remove low MAF loci
      if (ncol(gt_group_missing) == 0) next
    }

    if(allele_count_min>0){
      # for each loci, determine if the main or alternate allele has fewer values than the allele_count_min -- then remove
      loci_min_ac <- get_minor_allele_counts(gt_group_missing)
      passing_ac_loci <- which(loci_min_ac >= allele_count_min) # get indices
      gt_group_missing <- gt_group_missing[, ..passing_ac_loci] # remove low allele count loci
      if (ncol(gt_group_missing) == 0) next
    }
    
    group_out_list <- list() # make empty list to store group level data in 
    
    for (s in 1:length(sites)) { # for each site
      site <- sites[s] # get site name
      
      if(length(which(site_variable_group ==  site))<=1){
        next
      }
      
      gt_site <- gt_group_missing[which(site_variable_group ==  site), ] # get samples gt data
      
      n_site <- n
      if(is.null(n)){
        n_site <- nrow(gt_site)
        if(n_site <= 3) next # hard minimum of 3 samples per group
      }
      
      # loci <- ncol(gt_site) 
      
      # what to do if there are fewer samples than number of individuals requested?
      if(n_site<=nrow(gt_site)){ # if there are more or equal number of samples requested, proceed
        ind_index <- sample(1:nrow(gt_site), n_site, replace = FALSE)
        gt_site_subset <- gt_site[ind_index, ]
        
      } else{ # if there are fewer than requested, skip site
        next
      }
      
      not_na_loci <- which(colSums(is.na(gt_site_subset))/nrow(gt_site_subset)<1)
      gt_site_subset2 <- gt_site_subset[,..not_na_loci]
      
      # # TRIAL
      # not_missing_loci <- which(colMeans(is.na(gt_site_subset2)) <=  max_missingness)  # determine which loci pass missingness threshold
      # gt_site_subset2 <- gt_site_subset2[, ..not_missing_loci] # filter to remove high missing loci
      # if (ncol(gt_site_subset2) == 0) next
      # 
      # if(maf>0){
      #   loci_mafs <- get_minor_allele_frequencies(gt_site_subset2) # get minor allele frequencies of each locus
      #   passing_maf_loci <- which(loci_mafs >= maf) # get indices
      #   gt_site_subset2 <- gt_site_subset2[, ..passing_maf_loci] # remove low MAF loci
      #   if (ncol(gt_site_subset2) == 0) next
      # }
      # 
      # if(allele_count_min>0){
      #   # for each loci, determine if the main or alternate allele has fewer values than the allele_count_min -- then remove
      #   loci_min_ac <- get_minor_allele_counts(gt_site_subset2)
      #   passing_ac_loci <- which(loci_min_ac >= allele_count_min) # get indices
      #   gt_site_subset2 <- gt_site_subset2[, ..passing_ac_loci] # remove low allele count loci
      #   if (ncol(gt_site_subset2) == 0) next
      # }

      td_out <- calculate_pi_theta_d(gt_site_subset2)
      
      site_results <- c(
        "group" = group,
        "site" = site,
        "loci" = ncol(gt_site_subset2),
        "n" = nrow(gt_site_subset2), 
        "pi" = td_out$pi %>% as.numeric(), 
        "S" = td_out$S %>% as.numeric(),
        "theta" = td_out$theta %>% as.numeric(),
        "D" = td_out$D %>% as.numeric(),
        "real_missingness" = sum(is.na(gt_site_subset2))/(ncol(gt_site_subset2)*nrow(gt_site_subset2))
      )
      
      group_out_list[[s]] <- site_results
    }
    
    # out_list[[group]] <- bind_rows(group_out_list)
    if (length(group_out_list) > 0) {
      out_list[[group]] <- bind_rows(group_out_list)
      }
    }
    
    # out_matrix <- bind_rows(out_list)
    # out_data <- out_matrix
    if(length(out_list) == 0) {
      warning("No data passed filters")
      return(NULL)
    }
    
    out_matrix <- bind_rows(out_list)
  
    print("Process complete!")
    
    return(out_matrix)
}

