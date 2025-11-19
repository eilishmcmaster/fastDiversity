#' Calculate Tajima's D as in Pixy 
#' 
#' @description
#' This function calculates Tajima's D, Watterson's theta, and pi following the methods in Pixy, which accounts for missing data to produce unbiased estimates. For further details see https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13326 and https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.14104 
#' 
#' @param gt Genotype matrix where individuals are rows and loci are columns, coded as 
#'   0=aa, 1=aA, 2=AA. Note these methods are designed to have fixed/invariant loci 
#'   included (https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13326). 
tajimas_d <- function(gt) {
  
  site_stats <- apply(gt, 2, function(col) {
    nonmiss <- !is.na(col)
    k <- sum(nonmiss)   # diploid individuals
    n_s <- 2 * k  # haploid sample size
    # note this PIXY version includes fixed sites (where all alleles are the same) which is not the case in Tajima 1989
    
    if (n_s < 2) {
      return(c(n_s = n_s,
               mpd = 0,
               segregating = 0))
    }
    
    x <- sum(col[nonmiss])  # alt allele count
    
    if (x == 0 || x == n_s) { # if there are no alt alleles or alt allele is fixed, return 0
      return(c(n_s = n_s,
               mpd = 0,
               segregating = 0))
    }
    
    num_diff_pairs <- x * (n_s - x) # number of different haploid pairs (total differences Korunes + Samuk 2021)
    num_pairs <- n_s * (n_s - 1) / 2 # total number of haploid pairs (no. comparisons Korunes + Samuk 2021)
    
    return(c(
      n_s = n_s,
      num_diff_pairs = num_diff_pairs, 
      num_pairs = num_pairs, 
      segregating = 1
    ))
  })
  
  site_stats <- as.data.frame(t(site_stats))
  
  pi <- sum(site_stats$num_diff_pairs, na.rm = TRUE)/ sum(site_stats$num_pairs, na.rm = TRUE) # pi raw described in Korunes + Samuk 2021 
  
  S <- sum(site_stats$segregating) # number of sites/loci that are segregating/polymorphic
  
  # wattersons theta used by pixy described in Bailey, Stevison, Samuk 2025
  # pixy groups segregating sites by n_s then W = sum(S_n / a1_n)
  
  theta_w <- 0
  segregating_sites <- site_stats[site_stats$segregating == 1, ]
  
  if (nrow(segregating_sites) > 0) { 
    # group by haploid sample size
    for (n_s in unique(segregating_sites$n_s)) { # group sites/loci that share the same number of haploid genotypes (handles missingness)
      S_n <- sum(segregating_sites$n_s == n_s)
      a1_n <- sum(1 / (1:(n_s - 1)))
      theta_w <- theta_w + S_n / a1_n
    }
  }
  
  # calcs for tajimas D
  n <- mean(site_stats$n_s, na.rm = TRUE) # mean haploid count (number of haploid sequences)
  if (n < 3){ # if there are less than 3 haploid sequences D is NA
    return(list(pi = pi, S = S, theta_w = theta_w, Tajimas_D = NA,
                per_site_stats = site_stats))
  }

  a1 <- sum(1 / (1:(n - 1))) # harmonic numbers from Tajima 1989
  a2 <- sum(1 / ((1:(n - 1)))^2)

  b1 <- (n + 1) / (3 * (n - 1))
  b2 <- 2 * (n^2 + n + 3) / (9 * n * (n - 1))

  c1 <- b1 - (1 / a1)
  c2 <- b2 - ((n + 2) / (a1 * n)) + (a2 / a1^2)

  e1 <- c1 / a1
  e2 <- c2 / (a1^2 + a2)

  V <- e1 * S + e2 * S * (S - 1) # variance
  if (V == 0) {
    D <- NA
  } else {
    D <- (pi - theta_w) / sqrt(V)
  }

  # ## original Tajimas D code and example:
  # # pi <- 15.38 # average pairwise nucleotide differences between individuals
  # # S <- 45 # number of segregating sites
  # # n <- 7 # number of individuals
  # 
  # a1 <- sum(1 / (1:(n - 1))) # should be 2.45
  # a2 <- sum(1 / ((1:(n - 1))^2)) # 1.49 
  # e1 <- (1/a1) * ((n+1)/(3*(n-1)) - 1/a1) #0.0148 
  # c <- ((2*(n^2 + n + 3))/(9*n*(n-1))) - ((n+2)/(a1*n)) + (a2/(a1^2)) #0.0358 
  # e2 <- c/(a1^2 + a2) #0.00478 
  # D <- (pi - S/a1)/ sqrt(e1*S + e2*S*(S-1)) #-0.938 
  ##

  return(list(
    pi = pi,
    S = S,
    theta_w = theta_w,
    Tajimas_D = D,
    per_site_stats = site_stats
  ))
}