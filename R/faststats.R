#' Calculate basic diversity statistics
#' 
#' This function calculates allelic richness, heterozygosity, and inbreeding coefficients 
#' for genetic groups within specified sites. Data is filtered based on genetic group and 
#' site IDs, and parameters such as minimum number of individuals per site, minimum number 
#' of loci, minimum minor allele frequency, and maximum missingness per locus.
#' 
#' @param gt Genotype matrix where individuals are rows and loci are columns, coded as 
#'           0=aa, 1=aA, 2=AA.
#' @param site_variable Vector of group IDs (sites) for each individual (same length as 
#'                      number of rows in gt).
#' @param genetic_group_variable Vector of genetic group IDs for each individual (same 
#'                               length as number of rows in gt). Can represent species, 
#'                               subspecies, or other biologically relevant structure groups.
#' @param minimum_n Minimum number of individuals required per site (sites with fewer 
#'                  individuals will be filtered out).
#' @param minimum_loci Minimum number of loci required to proceed with analysis for each 
#'                    genetic group (genetic groups with fewer loci will be skipped).
#' @param maf Minimum minor allele frequency filter.
#' @param max_missingness Maximum proportion of missing data allowed per locus (loci with 
#'                        higher missingness will be filtered out).
#' @return A data frame containing diversity statistics (allelic richness, heterozygosity, 
#'         inbreeding coefficients) per site per genetic group.
#' @examples
#' data(example_data)
#' # Calculate observed heterozygosity (Ho), expected heterozygosity (He), unbiased 
#' # expected heterozygosity (uHe), inbreeding coefficient (Fis), mean allelic richness (Ar),
#' # and rarefied allelic richness (rAr):
#' basicstats <- faststats(example_gt, example_meta$population,
#'                         example_meta$site, minimum_n=3, 
#'                         minimum_loci=50, maf=0.05, max_missingness=0.3)
#' @export
faststats <- function (gt, genetic_group_variable, site_variable, minimum_n = 3, 
                       minimum_loci = 50, maf = 0.05, max_missingness = 0.3, 
                       fis_ci = FALSE, boots = NULL, resample_n = NULL, fis_alpha=0.05) 
{
  gt <- data.table::as.data.table(gt, keep.rownames = FALSE)
  out_list <- list()
  
  genetic_group_freq <- table(genetic_group_variable)
  genetic_groups <- names(which(genetic_group_freq >= minimum_n))
  
  for (group in genetic_groups) {
    gt_group <- gt[which(genetic_group_variable == group), ]
    site_variable_group <- site_variable[which(genetic_group_variable ==  group)]
    not_missing_loci <- which(colMeans(is.na(gt_group)) <=  max_missingness) 
    
    gt_group_missing <- gt_group[, ..not_missing_loci]
    
    loci_mafs <- get_minor_allele_frequencies(gt_group_missing)
    passing_maf_loci <- which(loci_mafs >= maf)
    gt_group_missing_maf <- gt_group_missing[, ..passing_maf_loci]
    
    if (length(passing_maf_loci) < minimum_loci) {
      print(paste(group, "does not have enough loci (",
                  length(passing_maf_loci), "loci: minimum is,",
                  minimum_loci, ")"))
      print("Proceeding to next genetic group")
      (next)()
    }
    site_freq <- table(site_variable_group)
    sites <- names(which(site_freq >= minimum_n))
    # group_out_matrix <- matrix(NA, length(sites), 10)
    group_out_list <- list()
    
    for (s in 1:length(sites)) {
      site <- sites[s]
      gt_site <- gt_group_missing_maf[which(site_variable_group == 
                                              site), ]
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
        "uFis" = uFis,
        "loci" = loci,
        "n" = n
      )
      
      #### bootstrapped Ho, He, Fis ####
      if(isTRUE(fis_ci) || !is.null(boots)){ # run if fis_ci is TRUE or value is supplied for boots
        boot_stats <- calculate_boot_stats(gt_site, boots, resample_n, fis_alpha=0.05)

        site_results <- c(site_results, boot_stats)
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
    print(paste(group, "complete"))
  }
  out_matrix <- bind_rows(out_list)
  print("Process complete!")
  return(out_matrix)
}
