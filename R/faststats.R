#' Calculate basic diversity statistics
#' 
#' @description
#' This function calculates allelic richness, heterozygosity, and inbreeding coefficients 
#' for specified sites within genetic groups. Data is filtered based on genetic group (minimum
#' number of loci, minimum minor allele frequency, and maximum missingness per locus) 
#' and minimum number of individuals per site.
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
#' @param minimum_loci Minimum number of loci required to proceed with analysis for each 
#'   genetic group (genetic groups with fewer loci will be skipped).
#' @param maf Minimum minor allele frequency filter.
#' @param max_missingness Maximum proportion of missing data allowed per locus (loci with 
#'   higher missingness will be filtered out).
#' @param get_CI_loci_resampling Get confidence intervals for Ho, He, Fis, using bootstrapping by resampling loci (as in hierfstat::boot.ppfis)
#' @param get_CI_individual_resampling Get confidence intervals for Ho, He, Fis, using bootstrapping by resampling individuals (as in diveRsity::basicstats)
#' @param boots Number of bootstraps to get confidence intervals.
#' @param resample_n Number of samples to resample per site if bootstrapping, defaults to site n.
#' @param CI_alpha Alpha value for confidence intervals, defaults to 0.05 (2.5\%, 97.5\% CIs).
#' @param run_HWE_test Specify whether to statistically test if allele frequencies differ from Hardy-Weinberg equilibrium.
#' @param return_locus_stats Specify whether to return U-score and LLR p values for each locus from HWE testing.
#' 
#' @return A data frame with the following columns:
#'   \item{group}{Genetic group ID (e.g., species, subspecies).}
#'   \item{site}{Site ID corresponding to the sample location.}
#'   \item{Ar}{Allelic richness per site per group, standardized for sample size.}
#'   \item{Ho}{Observed heterozygosity.}
#'   \item{He}{Expected heterozygosity under Hardy-Weinberg equilibrium within the site_variable.}
#'   \item{uHe}{Unbiased expected heterozygosity (corrected for small sample size).}
#'   \item{Fis}{Inbreeding coefficient (Fis = 1 - Ho / He).}
#'   \item{uF}{Unbiased inbreeding coefficient (Fis = 1 - Ho / uHe).}
#'   \item{loci}{Number of loci used in the analysis for that group-site combination.}
#'   \item{n}{Number of individuals sampled at the site for the given group.}
#'   \item{*_ci.*\%, *_ci.*\%}{Lower and upper confidence intervals for Fis, Ho, and He, generated via bootstrapping. Column prefixes indicate the statistic (e.g., \code{loci_f.2.5\%}, \code{loci_he.97.5\%}).}
#'   \item{LLR_fisher_p}{P-value from Fisher's method combining Hardy-Weinberg equilibrium log-likelihood ratio (LLR) tests across loci.}
#'   \item{mean_U_score}{Mean U-score across loci from Hardy-Weinberg tests. Negative values indicate heterozygosity excess; positive values indicate homozygosity excess.}
#' 
#' @import data.table HWxtest dplyr stats
#' @seealso [resampling_individuals_CI(), resampling_loci_CI()]
#' 
#' 
#' @export
faststats <- function (gt, genetic_group_variable, site_variable, minimum_n = 3, 
                       minimum_loci = 50, maf = 0.05, max_missingness = 0.3, 
                       get_CI_loci_resampling = FALSE, loci_resample_n = NULL, 
                       get_CI_individual_resampling = FALSE, individual_resample_n = NULL, 
                       boots = 100, CI_alpha=0.05,
                       run_HWE_test = FALSE, return_locus_stats = FALSE) 
{
  locus_names <- colnames(gt)
  gt <- data.table::as.data.table(gt, keep.rownames = FALSE) # make genotype data into data.table
  setnames(gt, locus_names)
  
  out_list <- list() # make empty list to put data in 
  
  if(isTRUE(return_locus_stats)){ # if locus stats should be returned
    locus_out_list <- list()
  }
  
  genetic_group_freq <- table(genetic_group_variable) # get number of samples per genetic group
  genetic_groups <- names(which(genetic_group_freq >= minimum_n)) # get the genetic groups with n>=minimum
  
  for (group in genetic_groups) { # for each genetic group
    gt_group <- gt[which(genetic_group_variable == group), ] # filter samples
    site_variable_group <- site_variable[which(genetic_group_variable ==  group)] # get vector of sites for samples in the group
    not_missing_loci <- which(colMeans(is.na(gt_group)) <=  max_missingness)  # determine which loci pass missingness threshold
    
    gt_group_missing <- gt_group[, ..not_missing_loci] # filter to remove high missing loci
    
    loci_mafs <- get_minor_allele_frequencies(gt_group_missing) # get minor allele frequencies of each locus
    passing_maf_loci <- which(loci_mafs >= maf) # get indices 
    gt_group_missing_maf <- gt_group_missing[, ..passing_maf_loci] # remove low MAF loci
    
    if (length(passing_maf_loci) < minimum_loci) {
      print(paste(group, "does not have enough loci (",
                  length(passing_maf_loci), "loci: minimum is,",
                  minimum_loci, ")"))
      print("Proceeding to next genetic group")
      (next)()
    }
    
    site_freq <- table(site_variable_group) # get number of samples per site
    sites <- names(which(site_freq >= minimum_n)) # remove sites less than threshold n

    group_out_list <- list() # make empty list to store group level data in 
    if(isTRUE(return_locus_stats)){ # if locus stats should be returned
      group_locus_out_df <- data.table()
    }
    
    for (s in 1:length(sites)) { # for each site
      site <- sites[s] # get site name
      gt_site <- gt_group_missing_maf[which(site_variable_group ==  site), ] # get samples gt data
      
      # calculate basic stats
      Ho <- calculate_Ho(gt_site) %>% round(.,3)
      Hes <- calculate_Hes(gt_site) %>% round(.,3)
      He <- mean(Hes, na.rm = TRUE) %>% round(.,3)
      ns <- colSums(!is.na(gt_site))
      uHe <- calculate_uHe(ns, Hes) %>% round(.,3)
      Fis <- (1 - (Ho/He) ) %>% round(.,3)
      uFis <- (1 - (Ho/uHe)) %>% round(.,3)
      loci <- ncol(gt_site) 
      Ar <- calculate_Ar(gt_site) %>% round(.,3)
      n <- nrow(gt_site)
      
      site_results <- c(
        "group" = group,
        "site" = site,
        "Ar" = Ar,
        "Ho" = Ho,
        "He" = He,
        "uHe" = uHe,
        "Fis" = Fis,
        "uF" = uFis,
        "loci" = loci,
        "n" = n
      )
      
      #### bootstrapped Ho, He, Fis ####
      if(isTRUE(get_CI_individual_resampling)){ # run if get_CI is TRUE
        boot_stats <- resampling_individuals_CI(gt_site, boots, individual_resample_n, CI_alpha=CI_alpha)
        site_results <- c(site_results, boot_stats)
      }
      
      if(isTRUE(get_CI_loci_resampling)){ # run if get_CI is TRUE
        boot_stats <- resampling_loci_CI(gt_site, boots, loci_resample_n, CI_alpha=CI_alpha)
        site_results <- c(site_results, boot_stats)
      }
      
      #### HWE test for loci ####
      if(isTRUE(run_HWE_test)){ # run HWE tests 
        HWE_test_loci <- HWEtest_of_loci(gt_site)
        site_results <- c(site_results, (HWE_test_loci[1:2]%>% unlist))
      }
      
      if(isTRUE(return_locus_stats)){ # if locus data is requested, add to the table
        group_locus_out_df <- bind_rows(group_locus_out_df, data.table(HWE_test_loci[[3]], 'site'=site))
      }
      
      group_out_list[[s]] <- site_results
    }
    
    # # rarefied AR
    # count_matrix <- t(sapply(split(seq_along(site_variable_group), 
    #                                site_variable_group), function(idx) colSums(gt_group_missing_maf[idx, 
    #                                ], na.rm = TRUE))) %>% as.matrix()
    # n_matrix <- t(sapply(split(seq_along(site_variable_group), 
    #                            site_variable_group), function(idx) colSums(!is.na(gt_group_missing_maf[idx, 
    #                            ])))) %>% as.matrix()
    # a_matrix <- 2 * n_matrix
    # min_n_per_locus <- apply(n_matrix, 2, min)
    # n <- min_n_per_locus * 2
    # count_matrix2 <- a_matrix - count_matrix
    # all_ar <- matrix(NA, nrow(n_matrix), ncol(n_matrix))
    # for (i in 1:ncol(n_matrix)) {
    #   for (j in 1:nrow(n_matrix)) {
    #     samp <- choose(a_matrix[j, i] - c(count_matrix[j, 
    #                                                    i], count_matrix2[j, i]), n[i])/choose(a_matrix[j, 
    #                                                                                                    i], n[i])
    #     samp[is.na(samp)] <- 0
    #     all_ar[j, i] <- sum(1 - samp)
    #   }
    # }
    # mn_ar <- rowMeans(all_ar, na.rm = TRUE) %>% round(., 
    #                                                   3)
    # sd_ar <- apply(all_ar, 1, sd, na.rm = TRUE) %>% round(., 
    #                                                       3)
    # group_out_matrix <- cbind(group_out_matrix, mn_ar, sd_ar)
    # out_matrix <- rbind(out_matrix, group_out_matrix)
    
    out_list[[group]] <- group_out_list
    
    if(isTRUE(return_locus_stats)){ # if locus data is requested, add to the table
      locus_out_list[[group]] <- group_locus_out_df 
    }
    print(paste(group, "complete"))
  }
  
  out_matrix <- bind_rows(out_list)
  
  if(isTRUE(return_locus_stats)){ # if locus data is requested, add to the table
    out_data <- list('summary'= out_matrix, 'locus_HWE' = locus_out_list)
  }else{
    out_data <- out_matrix
    }
  
  print("Process complete!")
  
  return(out_data)
  
}
