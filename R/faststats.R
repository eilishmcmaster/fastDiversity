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
#' @param loci_resample_n Number of loci to resample if bootstrapping, defaults to site n.
#' @param individual_resample_n Number of samples to resample if bootstrapping, defaults to site n.
#' @param CI_alpha Alpha value for confidence intervals, defaults to 0.05 (2.5\%, 97.5\% CIs).
#' @param run_HWE_test Specify whether to statistically test if allele frequencies differ from Hardy-Weinberg equilibrium.
#' @param return_locus_stats Specify whether to return p values for each locus from HWE testing.
#'   Setting this to TRUE also runs the HWE test (\code{run_HWE_test = TRUE}).
#' 
#' @return A data frame with the following columns (or, if \code{return_locus_stats = TRUE},
#'   a list with this data frame as \code{summary} and per-locus HWE p-values as \code{locus_HWE}):
#'   \item{group}{Genetic group ID (e.g., species, subspecies).}
#'   \item{site}{Site ID corresponding to the sample location.}
#'   \item{Ar}{Rarefied allelic richness: the expected number of alleles per locus in a sample of
#'     the smallest size observed at that locus among the sites retained for the genetic group,
#'     averaged across loci (see \code{\link{calculate_rarified_Ar}}). Comparable among sites within
#'     a genetic group.}
#'   \item{Ho}{Observed heterozygosity.}
#'   \item{He}{Expected heterozygosity under Hardy-Weinberg equilibrium within the site_variable.}
#'   \item{uHe}{Unbiased expected heterozygosity (corrected for small sample size).}
#'   \item{Fis}{Inbreeding coefficient (Fis = 1 - Ho / He).}
#'   \item{Gst}{Site-level gene diversity differentiation in the sense of Nei (1973),
#'     (Ht - Hs) / Ht, where Hs is the site's He and Ht is the He of all individuals in the
#'     genetic group.}
#'   \item{Fst}{Population-specific Fst of Weir & Goudet (2017) for the site, relative to the other
#'     sites retained for the genetic group (see \code{\link{calculate_Fst}}). NA if the group
#'     has only one retained site.}
#'   \item{uFis}{Unbiased inbreeding coefficient (uFis = 1 - Ho / uHe).}
#'   \item{loci}{Number of loci used in the analysis for that group-site combination.}
#'   \item{n}{Number of individuals sampled at the site for the given group.}
#'   \item{*_ci.*\%, *_ci.*\%}{Lower and upper confidence intervals for Fis, Ho, and He, generated via bootstrapping. Column prefixes indicate the statistic (e.g., \code{loci_f.2.5\%}, \code{loci_he.97.5\%}).}
#'   \item{HWE_p}{P-value from Fisher's method combining Hardy-Weinberg equilibrium tests across loci (HardyWeinberg::HWExact).}
#' 
#' @references
#' Nei, M. (1973). Analysis of gene diversity in subdivided populations. Proceedings of the
#' National Academy of Sciences, 70(12), 3321--3323. https://doi.org/10.1073/pnas.70.12.3321
#'
#' Weir, B. S., & Goudet, J. (2017). A unified characterization of population structure and
#' relatedness. Genetics, 206(4), 2085--2103. https://doi.org/10.1534/genetics.116.198424
#'
#' Kalinowski, S. T. (2004). Counting alleles with rarefaction: private alleles and hierarchical
#' sampling designs. Conservation Genetics, 5(4), 539--543.
#' https://doi.org/10.1023/B:COGE.0000041021.91777.1a
#'
#' @seealso [resampling_individuals_CI(), resampling_loci_CI(), calculate_Fst(),
#'   calculate_rarified_Ar()]
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
    if(isTRUE(return_locus_stats) && !isTRUE(run_HWE_test)){ # locus stats come from the HWE test
      message("return_locus_stats = TRUE requires the HWE test; setting run_HWE_test = TRUE.")
      run_HWE_test <- TRUE
    }
    
    locus_names <- colnames(gt)
    gt <- data.table::as.data.table(gt, keep.rownames = FALSE) # make genotype data into data.table
    data.table::setnames(gt, locus_names)
    
    out_list <- list() # make empty list to put data in 
    
    if(isTRUE(return_locus_stats)){ # if locus stats should be returned
      locus_out_list <- list()
    }
    
    genetic_group_freq <- table(genetic_group_variable) # get number of samples per genetic group
    genetic_groups <- names(which(genetic_group_freq >= minimum_n)) # get the genetic groups with n>=minimum
    
    for (group in genetic_groups) { # for each genetic group
      gt_group <- gt[which(genetic_group_variable == group), ] # filter samples
      site_variable_group <- site_variable[which(genetic_group_variable ==  group)] # get vector of sites for samples in the group
      #
      site_freq <- table(site_variable_group) # get number of samples per site
      sites <- names(which(site_freq >= minimum_n)) # remove sites less than threshold n
      
      if (length(sites) == 0) {
        message(sprintf("Skipping group '%s': no sites meet n >= %d.", 
                        group, minimum_n)) 
        next  }
      #
      not_missing_loci <- which(colMeans(is.na(gt_group)) <=  max_missingness)  # determine which loci pass missingness threshold
      
      if (length(not_missing_loci) == 0) {
        message(group, " has no loci with missingness <= ", max_missingness,
                ". Proceeding to next genetic group")
        next
      }
      gt_group_missing <- gt_group[, ..not_missing_loci] # filter to remove high missing loci
      
      loci_mafs <- get_minor_allele_frequencies(gt_group_missing) # get minor allele frequencies of each locus
      passing_maf_loci <- which(loci_mafs >= maf) # get indices 
      gt_group_missing_maf <- gt_group_missing[, ..passing_maf_loci] # remove low MAF loci
      
      if (length(passing_maf_loci) < max(minimum_loci, 1)) { # at least one locus is always needed
        message(group, " does not have enough loci (", length(passing_maf_loci),
                " loci: minimum is ", minimum_loci, "). Proceeding to next genetic group")
        next
      }
      
      Ht <- mean(calculate_Hes(gt_group_missing_maf), na.rm = TRUE) # expected heterozygosity across all individuals and averaged across loci
      # rarefied allelic richness and Weir & Goudet Fst are relative to the other retained sites
      retained <- which(site_variable_group %in% sites)
      gt_retained <- as.matrix(gt_group_missing_maf[retained, ])
      sites_retained <- as.character(site_variable_group[retained])
      Ar_sites <- rowMeans(calculate_rarified_Ar(gt_retained, sites_retained), na.rm = TRUE)
      if (length(sites) > 1) {
        Fst_sites <- calculate_Fst(gt_retained, sites_retained)
      } else {
        Fst_sites <- setNames(NA_real_, sites)
      }
      
      group_out_list <- list() # make empty list to store group level data in 
      if(isTRUE(return_locus_stats)){ # if locus stats should be returned
        group_locus_out_df <- data.table()
      }
      
      for (s in 1:length(sites)) { # for each site
        site <- sites[s] # get site name
        gt_site <- gt_group_missing_maf[which(site_variable_group ==  site), ] # get samples gt data
        
        # calculate basic stats (at full precision; rounded to 4 decimal places in the output)
        Ho <- calculate_Ho(gt_site)
        Hes <- calculate_Hes(gt_site)
        He <- mean(Hes, na.rm = TRUE)
        Gst <- (Ht-He)/(Ht) # Nei 1973 gene diversity differentiation for this site
        Fst <- Fst_sites[[site]] # Weir & Goudet 2017 population-specific Fst
        ns <- colSums(!is.na(gt_site))
        uHe <- calculate_uHe(ns, Hes)
        Fis <- 1 - (Ho/He)
        uFis <- 1 - (Ho/uHe)
        loci <- ncol(gt_site) 
        Ar <- Ar_sites[[site]]
        n <- nrow(gt_site)
        
        site_results <- c(
          "group" = group,
          "site" = site,
          "Ar" = round(Ar, 4),
          "Ho" = round(Ho, 4),
          "He" = round(He, 4),
          "uHe" = round(uHe, 4),
          "Fis" = round(Fis, 4),
          "Gst" = round(Gst, 4),
          "Fst" = round(Fst, 4),
          "uFis" = round(uFis, 4),
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
          HWE_test_loci <- MC_HWEtest_of_loci(gt_site)
          HWE_overall_p <- HWE_test_loci[[1]]
          names(HWE_overall_p) <- 'HWE_p'
          site_results <- c(site_results, HWE_overall_p)
        }
        
        if(isTRUE(return_locus_stats)){ # if locus data is requested, add to the table
          group_locus_out_df <- bind_rows(group_locus_out_df, c('site'=site, HWE_test_loci[[2]]))
        }
        
        group_out_list[[s]] <- site_results
      }
      
      out_list[[group]] <- bind_rows(group_out_list)
      
      if(isTRUE(return_locus_stats)){ # if locus data is requested, add to the table
        locus_out_list[[group]] <- group_locus_out_df 
      }
      message(group, " complete")
    }
    
    out_matrix <- bind_rows(out_list)
    
    # site_results is built with c(), which makes every value character; restore numeric columns
    numeric_cols <- setdiff(names(out_matrix), c("group", "site"))
    out_matrix[numeric_cols] <- lapply(out_matrix[numeric_cols], as.numeric)
    
    if(isTRUE(return_locus_stats)){ # if locus data is requested, add to the table
      out_data <- list('summary'= out_matrix, 'locus_HWE' = locus_out_list)
    }else{
      out_data <- out_matrix
    }
    
    message("Process complete!")
    
    return(out_data)
    
  }
