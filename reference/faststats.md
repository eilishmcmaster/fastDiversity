# Calculate basic diversity statistics

This function calculates allelic richness, heterozygosity, and
inbreeding coefficients for specified sites within genetic groups. Data
is filtered based on genetic group (minimum number of loci, minimum
minor allele frequency, and maximum missingness per locus) and minimum
number of individuals per site.

## Usage

``` r
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
```

## Arguments

- gt:

  Genotype matrix where individuals are rows and loci are columns, coded
  as 0=aa, 1=aA, 2=AA.

- genetic_group_variable:

  Vector of genetic group IDs for each individual (same length as number
  of rows in gt). Can represent species, subspecies, or other
  biologically relevant structure groups.

- site_variable:

  Vector of group IDs (sites) for each individual (same length as number
  of rows in gt).

- minimum_n:

  Minimum number of individuals required per site (sites with fewer
  individuals will be filtered out).

- minimum_loci:

  Minimum number of loci required to proceed with analysis for each
  genetic group (genetic groups with fewer loci will be skipped).

- maf:

  Minimum minor allele frequency filter.

- max_missingness:

  Maximum proportion of missing data allowed per locus (loci with higher
  missingness will be filtered out).

- get_CI_loci_resampling:

  Get confidence intervals for Ho, He, Fis, using bootstrapping by
  resampling loci (as in hierfstat::boot.ppfis)

- loci_resample_n:

  Number of loci to resample if bootstrapping, defaults to site n.

- get_CI_individual_resampling:

  Get confidence intervals for Ho, He, Fis, using bootstrapping by
  resampling individuals (as in diveRsity::basicstats)

- individual_resample_n:

  Number of samples to resample if bootstrapping, defaults to site n.

- boots:

  Number of bootstraps to get confidence intervals.

- CI_alpha:

  Alpha value for confidence intervals, defaults to 0.05 (2.5%, 97.5%
  CIs).

- run_HWE_test:

  Specify whether to statistically test if allele frequencies differ
  from Hardy-Weinberg equilibrium.

- return_locus_stats:

  Specify whether to return p values for each locus from HWE testing.

## Value

A data frame with the following columns:

- group:

  Genetic group ID (e.g., species, subspecies).

- site:

  Site ID corresponding to the sample location.

- Ar:

  Allelic richness per site per group, standardized for sample size.

- Ho:

  Observed heterozygosity.

- He:

  Expected heterozygosity under Hardy-Weinberg equilibrium within the
  site_variable.

- uHe:

  Unbiased expected heterozygosity (corrected for small sample size).

- Fis:

  Inbreeding coefficient (Fis = 1 - Ho / He).

- Fst:

  Fst/Gst Inbreeding due to population subdivision (Fst = (Ht-Hs)/Ht
  \[W+C 1984\])

- uF:

  Unbiased inbreeding coefficient (Fis = 1 - Ho / uHe).

- loci:

  Number of loci used in the analysis for that group-site combination.

- n:

  Number of individuals sampled at the site for the given group.

- \*\_ci.\*%, \*\_ci.\*%:

  Lower and upper confidence intervals for Fis, Ho, and He, generated
  via bootstrapping. Column prefixes indicate the statistic (e.g.,
  `loci_f.2.5%`, `loci_he.97.5%`).

- HW_fisher_p:

  P-value from Fisher's method combining Hardy-Weinberg equilibrium
  tests across loci (HardyWeinberg::HWExact).

## See also

\[resampling_individuals_CI(), resampling_loci_CI()\]
