#' Calculate Tajima's D, pi, and Watterson's theta
#'
#' @description
#' Calculates Tajima's D, nucleotide diversity (pi), and Watterson's theta
#' from a biallelic genotype matrix, following the approach of Korunes & Samuk
#' (2021) to account for missing data. Non-segregating and monomorphic sites
#' are removed prior to calculation.
#'
#' @param gt Genotype matrix where individuals are rows and loci are columns,
#'   coded as 0=aa, 1=aA, 2=AA.
#'
#' @return A list with elements:
#'   \item{tajima_d}{Tajima's D statistic}
#'   \item{num_sites}{Number of sites used after filtering}
#'   \item{raw_pi}{Nucleotide diversity (pi)}
#'   \item{watterson_theta}{Watterson's theta}
#'   \item{d_stdev}{Standard deviation of the denominator}
#'
#' @references
#' Korunes, K. L., & Samuk, K. (2021). pixy: Unbiased estimation of
#' nucleotide diversity and divergence in the presence of missing data.
#' Molecular Ecology Resources, 21(4), 1359–1368.
#' https://doi.org/10.1111/1755-0998.13326
#'
#' Tajima, F. (1989). Statistical method for testing the neutral mutation
#' hypothesis by DNA polymorphism. Genetics, 123(3), 585–595.
#' https://doi.org/10.1093/genetics/123.3.585
#'
#' @export
calculate_pi_theta_d <- function(gt) {
  snp_mat <- as.matrix(gt)
  
  # remove loci that are not segregating (all 0 or all 2) or are empty
  obs_counts <- colSums(!is.na(snp_mat))
  all_ref <- colSums(snp_mat == 0, na.rm = TRUE) == obs_counts
  all_alt <- colSums(snp_mat == 2, na.rm = TRUE) == obs_counts
  keep <- obs_counts > 0 & !(all_ref | all_alt)
  snp_mat <- snp_mat[, keep]
  
  # allele counts per locus
  ref_counts   <- colSums(2 - snp_mat, na.rm = TRUE)
  alt_counts   <- colSums(snp_mat, na.rm = TRUE)
  total_counts <- ref_counts + alt_counts
  
  valid        <- total_counts > 1
  alt_counts   <- alt_counts[valid]
  ref_counts   <- ref_counts[valid]
  total_counts <- total_counts[valid]
  
  num_sites <- sum(total_counts > 0)
  
  # pi
  p         <- alt_counts / total_counts
  q         <- ref_counts / total_counts
  pi_values <- 2 * p * q * total_counts / (total_counts - 1)
  pi_values[total_counts <= 1] <- 0
  raw_pi    <- sum(pi_values, na.rm = TRUE)
  
  # Watterson's theta
  variant_idx   <- alt_counts > 0 & alt_counts < total_counts
  variant_total <- total_counts[variant_idx]
  a1_vec        <- 1 / (seq_len(max(variant_total) - 1))
  watterson_theta <- sum(1 / sapply(variant_total,
                                    function(n) sum(a1_vec[1:(n - 1)])))
  
  # Tajima's D
  n  <- round(mean(total_counts))
  s  <- length(variant_total)
  a1 <- sum(1 / seq(1, n - 1))
  a2 <- sum(1 / (seq(1, n - 1)^2))
  b1 <- (n + 1) / (3 * (n - 1))
  b2 <- 2 * (n^2 + n + 3) / (9 * n * (n - 1))
  c1 <- b1 - (1 / a1)
  c2 <- b2 - ((n + 2) / (a1 * n)) + (a2 / (a1^2))
  e1 <- c1 / a1
  e2 <- c2 / (a1^2 + a2)
  d_stdev  <- sqrt((e1 * s) + (e2 * s * (s - 1)))
  tajima_d <- if (d_stdev > 0) (raw_pi - watterson_theta) / d_stdev else NA
  
  return(list(
    tajima_d        = tajima_d,
    num_sites       = num_sites,
    raw_pi          = raw_pi,
    watterson_theta = watterson_theta,
    d_stdev         = d_stdev
  ))
}
