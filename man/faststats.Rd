% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/faststats.R
\name{faststats}
\alias{faststats}
\title{Calculate basic diversity statistics}
\usage{
faststats(
  gt,
  genetic_group_variable,
  site_variable,
  minimum_n = 3,
  minimum_loci = 50,
  maf = 0.05,
  max_missingness = 0.3,
  get_CI_loci_resampling = FALSE,
  loci_resample_n = NULL,
  get_CI_individual_resampling = FALSE,
  individual_resample_n = NULL,
  boots = 100,
  CI_alpha = 0.05,
  run_HWE_test = FALSE,
  return_locus_stats = FALSE
)
}
\arguments{
\item{gt}{Genotype matrix where individuals are rows and loci are columns, coded as 
0=aa, 1=aA, 2=AA.}

\item{genetic_group_variable}{Vector of genetic group IDs for each individual (same 
length as number of rows in gt). Can represent species, subspecies, or other 
biologically relevant structure groups.}

\item{site_variable}{Vector of group IDs (sites) for each individual (same length as 
number of rows in gt).}

\item{minimum_n}{Minimum number of individuals required per site (sites with fewer 
individuals will be filtered out).}

\item{minimum_loci}{Minimum number of loci required to proceed with analysis for each 
genetic group (genetic groups with fewer loci will be skipped).}

\item{maf}{Minimum minor allele frequency filter.}

\item{max_missingness}{Maximum proportion of missing data allowed per locus (loci with 
higher missingness will be filtered out).}

\item{get_CI_loci_resampling}{Get confidence intervals for Ho, He, Fis, using bootstrapping by resampling loci (as in hierfstat::boot.ppfis)}

\item{get_CI_individual_resampling}{Get confidence intervals for Ho, He, Fis, using bootstrapping by resampling individuals (as in diveRsity::basicstats)}

\item{boots}{Number of bootstraps to get confidence intervals.}

\item{CI_alpha}{Alpha value for confidence intervals, defaults to 0.05 (2.5\%, 97.5\% CIs).}

\item{run_HWE_test}{Specify whether to statistically test if allele frequencies differ from Hardy-Weinberg equilibrium.}

\item{return_locus_stats}{Specify whether to return U-score and LLR p values for each locus from HWE testing.}

\item{resample_n}{Number of samples to resample per site if bootstrapping, defaults to site n.}
}
\value{
A data frame with the following columns:
  \item{group}{Genetic group ID (e.g., species, subspecies).}
  \item{site}{Site ID corresponding to the sample location.}
  \item{Ar}{Allelic richness per site per group, standardized for sample size.}
  \item{Ho}{Observed heterozygosity.}
  \item{He}{Expected heterozygosity under Hardy-Weinberg equilibrium within the site_variable.}
  \item{uHe}{Unbiased expected heterozygosity (corrected for small sample size).}
  \item{Fis}{Inbreeding coefficient (Fis = 1 - Ho / He).}
  \item{uF}{Unbiased inbreeding coefficient (Fis = 1 - Ho / uHe).}
  \item{loci}{Number of loci used in the analysis for that group-site combination.}
  \item{n}{Number of individuals sampled at the site for the given group.}
  \item{*_ci.*\%, *_ci.*\%}{Lower and upper confidence intervals for Fis, Ho, and He, generated via bootstrapping. Column prefixes indicate the statistic (e.g., \code{loci_f.2.5\%}, \code{loci_he.97.5\%}).}
  \item{LLR_fisher_p}{P-value from Fisher's method combining Hardy-Weinberg equilibrium log-likelihood ratio (LLR) tests across loci.}
  \item{mean_U_score}{Mean U-score across loci from Hardy-Weinberg tests. Negative values indicate heterozygosity excess; positive values indicate homozygosity excess.}
}
\description{
This function calculates allelic richness, heterozygosity, and inbreeding coefficients 
for specified sites within genetic groups. Data is filtered based on genetic group (minimum
number of loci, minimum minor allele frequency, and maximum missingness per locus) 
and minimum number of individuals per site.
}
\seealso{
[resampling_individuals_CI(), resampling_loci_CI()]
}
