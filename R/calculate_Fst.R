#' Calculate population-specific Fst (Weir & Goudet 2017)
#'
#' This function calculates population-specific Fst for each site using the
#' allele-matching estimator of Weir & Goudet (2017). For each site, Fst
#' compares the probability that two alleles drawn from different individuals
#' within that site are identical (within-site matching) with the average
#' probability that two alleles drawn from different sites are identical
#' (between-site matching): \eqn{F_{ST,i} = (M_i - M_B) / (1 - M_B)}.
#'
#' This is the population-specific counterpart of Weir & Cockerham's (1984)
#' \eqn{\theta}, and gives the same values as the \code{Fs["Fst", ]} row of
#' \code{hierfstat::fs.dosage()}. Missing genotypes are handled as in hierfstat,
#' by averaging matching over the loci genotyped in both individuals of a pair.
#' Values are relative to the set of sites supplied, so they depend on which
#' sites are included.
#'
#' @param gt Genotype matrix where individuals are rows and loci are columns, coded as
#'           0=aa, 1=aA, 2=AA.
#' @param site_variable Vector of site IDs for each individual (same length as
#'                      number of rows in gt). At least two sites are required.
#' @return Named vector of Fst values, one per site.
#' @references
#' Weir, B. S., & Goudet, J. (2017). A unified characterization of population
#' structure and relatedness. Genetics, 206(4), 2085--2103.
#' https://doi.org/10.1534/genetics.116.198424
#'
#' Weir, B. S., & Cockerham, C. C. (1984). Estimating F-statistics for the
#' analysis of population structure. Evolution, 38(6), 1358--1370.
#' https://doi.org/10.2307/2408641
#' @examples
#' data(example_data)
#' keep <- example_meta$population == 1
#' calculate_Fst(example_gt[keep, ], example_meta$site[keep])
#' @export
calculate_Fst <- function(gt, site_variable){
  gt <- as.matrix(gt)
  storage.mode(gt) <- "double"
  site_variable <- as.character(site_variable)
  sites <- unique(site_variable)
  if (length(sites) < 2) {
    stop("calculate_Fst() requires at least two sites.")
  }

  # allele matching between each pair of individuals, averaged over loci typed in both
  # (Weir & Goudet 2017; same calculation as hierfstat::matching())
  not_na <- !is.na(gt)
  centred <- gt - 1
  centred[!not_na] <- 0
  Mij <- 0.5 * (1 + tcrossprod(centred) / tcrossprod(not_na))
  diag(Mij) <- NA # only pairs of different individuals

  site_idx <- lapply(sites, function(s) which(site_variable == s))

  # within-site matching, and mean matching between each pair of sites
  M_within <- vapply(site_idx, function(i) mean(Mij[i, i], na.rm = TRUE), numeric(1))
  site_pairs <- utils::combn(length(sites), 2)
  M_between <- mean(apply(site_pairs, 2, function(p) {
    mean(Mij[site_idx[[p[1]]], site_idx[[p[2]]]], na.rm = TRUE)
  }))

  Fst <- (M_within - M_between) / (1 - M_between)
  names(Fst) <- sites
  return(Fst)
}
