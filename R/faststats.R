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
#' # Calculate observed heterozygosity (Ho), expected heterozygosity (He), unbiased 
#' # expected heterozygosity (uHe), inbreeding coefficient (Fis), mean allelic richness (Ar),
#' # and rarefied allelic richness (rAr):
#' basicstats <- faststats(example_gt, example_meta$population,
#'                         example_meta$site, minimum_n=3, 
#'                         minimum_loci=50, maf=0.05, max_missingness=0.3)
#' @export
faststats <- function (gt, genetic_group_variable, site_variable, minimum_n = 3, 
                       minimum_loci = 50, maf = 0.05, max_missingness = 0.3) 
{
  out_matrix <- matrix(NA, 0, 12)
  colnames(out_matrix) <- c("genetic_group", "site", "Ar", 
                            "Ho", "He", "uHe", "Fis", "uFis", "loci", "n", "rAr", 
                            "sd_rAr")
  
  genetic_group_freq <- table(genetic_group_variable)
  genetic_groups <- names(which(genetic_group_freq >= minimum_n))
  
  for (group in genetic_groups) {
    gt_group <- gt[which(genetic_group_variable == group), 
    ]
    site_variable_group <- site_variable[which(genetic_group_variable == 
                                                 group)]
    not_missing_loci <- which(colMeans(is.na(gt_group)) <= 
                                max_missingness)
    gt_group_missing <- gt_group[, not_missing_loci]
    loci_mafs <- get_minor_allele_frequencies(gt_group_missing)
    passing_maf_loci <- which(loci_mafs >= maf)
    gt_group_missing_maf <- gt_group_missing[, passing_maf_loci]
    if (ncol(gt_group_missing_maf) < minimum_loci) {
      print(paste(group), "does not have enough loci (", 
            ncol(gt_group_missing_maf), "loci: minimum is,", 
            minimum_loci, ")")
      print("Proceeding to next genetic group")
      (next)()
    }
    site_freq <- table(site_variable_group)
    sites <- names(which(site_freq >= minimum_n))
    site_variable_group <- site_variable_group[which(site_variable_group %in% 
                                                       sites)]
    group_out_matrix <- matrix(NA, length(sites), 10)
    for (s in 1:length(sites)) {
      site <- sites[s]
      gt_site <- gt_group_missing_maf[which(site_variable_group == 
                                              site), ]
      Ho <- calculate_Ho(gt_site)
      Hes <- calculate_Hes(gt_site)
      He <- mean(Hes, na.rm = TRUE)
      ns <- colSums(!is.na(gt_site))
      uHe <- calculate_uHe(ns, Hes)
      Fis <- 1 - (Ho/He)
      uFis <- 1 - (Ho/uHe)
      loci <- ncol(gt_site)
      Ar <- calculate_Ar(gt_site)
      n <- nrow(gt_site)
      group_out_matrix[s, ] <- c(group, site, round(Ar, 
                                                    3), round(Ho, 3), round(He, 3), round(uHe, 3), 
                                 round(Fis, 3), round(uFis, 3), loci, n)
    }
    count_matrix <- t(sapply(split(seq_along(site_variable_group), 
                                   site_variable_group), function(idx) colSums(gt_group_missing_maf[idx, 
                                   ], na.rm = TRUE))) %>% as.matrix()
    n_matrix <- t(sapply(split(seq_along(site_variable_group), 
                               site_variable_group), function(idx) colSums(!is.na(gt_group_missing_maf[idx, 
                               ])))) %>% as.matrix()
    a_matrix <- 2 * n_matrix
    min_n_per_locus <- apply(n_matrix, 2, min)
    n <- min_n_per_locus * 2
    count_matrix2 <- a_matrix - count_matrix
    all_ar <- matrix(NA, nrow(n_matrix), ncol(n_matrix))
    for (i in 1:ncol(n_matrix)) {
      for (j in 1:nrow(n_matrix)) {
        samp <- choose(a_matrix[j, i] - c(count_matrix[j, 
                                                       i], count_matrix2[j, i]), n[i])/choose(a_matrix[j, 
                                                                                                       i], n[i])
        samp[is.na(samp)] <- 0
        all_ar[j, i] <- sum(1 - samp)
      }
    }
    mn_ar <- rowMeans(all_ar, na.rm = TRUE) %>% round(., 
                                                      3)
    sd_ar <- apply(all_ar, 1, sd, na.rm = TRUE) %>% round(., 
                                                          3)
    group_out_matrix <- cbind(group_out_matrix, mn_ar, sd_ar)
    out_matrix <- rbind(out_matrix, group_out_matrix)
    print(paste(group, "complete"))
  }
  print("Process complete!")
  return(out_matrix %>% as.data.frame())
}
