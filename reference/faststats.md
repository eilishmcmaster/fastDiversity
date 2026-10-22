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
  Setting this to TRUE also runs the HWE test (`run_HWE_test = TRUE`).

## Value

A data frame with the following columns (or, if
`return_locus_stats = TRUE`, a list with this data frame as `summary`
and per-locus HWE p-values as `locus_HWE`):

- group:

  Genetic group ID (e.g., species, subspecies).

- site:

  Site ID corresponding to the sample location.

- Ar:

  Rarefied allelic richness: the expected number of alleles per locus in
  a sample of the smallest size observed at that locus among the sites
  retained for the genetic group, averaged across loci (see
  [`calculate_rarified_Ar`](https://eilishmcmaster.github.io/fastDiversity/reference/calculate_rarified_Ar.md)).
  Comparable among sites within a genetic group.

- Ho:

  Observed heterozygosity.

- He:

  Expected heterozygosity under Hardy-Weinberg equilibrium within the
  site_variable.

- uHe:

  Unbiased expected heterozygosity (corrected for small sample size).

- Fis:

  Inbreeding coefficient (Fis = 1 - Ho / He).

- Gst:

  Site-level gene diversity differentiation in the sense of Nei (1973),
  (Ht - Hs) / Ht, where Hs is the site's He and Ht is the He of all
  individuals in the genetic group.

- Fst:

  Population-specific Fst of Weir & Goudet (2017) for the site, relative
  to the other sites retained for the genetic group (see
  [`calculate_Fst`](https://eilishmcmaster.github.io/fastDiversity/reference/calculate_Fst.md)).
  NA if the group has only one retained site.

- uFis:

  Unbiased inbreeding coefficient (uFis = 1 - Ho / uHe).

- loci:

  Number of loci used in the analysis for that group-site combination.

- n:

  Number of individuals sampled at the site for the given group.

- \*\_ci.\*%, \*\_ci.\*%:

  Lower and upper confidence intervals for Fis, Ho, and He, generated
  via bootstrapping. Column prefixes indicate the statistic (e.g.,
  `loci_f.2.5%`, `loci_he.97.5%`).

- HWE_p:

  P-value from Fisher's method combining Hardy-Weinberg equilibrium
  tests across loci (HardyWeinberg::HWExact).

## Details

**Monomorphic sites and comparability.** Statistics are calculated from
the loci that pass the filters within each genetic group. With `maf > 0`
(the default) these are polymorphic loci only, and most SNP-calling
pipelines output only polymorphic sites in any case. Ho, He, uHe and Ar
are therefore SNP-level values: they are comparable among sites within a
genetic group, which share the same loci, but they are not genome-wide
heterozygosity and should not be compared with values from other
datasets or from other genetic groups, which retain different loci. Fis,
uFis, Gst and Fst are ratios that do not change when monomorphic loci
are added (exactly for complete data, approximately with missing data),
so they can be compared more broadly.

## References

Nei, M. (1973). Analysis of gene diversity in subdivided populations.
Proceedings of the National Academy of Sciences, 70(12), 3321–3323.
https://doi.org/10.1073/pnas.70.12.3321

Weir, B. S., & Goudet, J. (2017). A unified characterization of
population structure and relatedness. Genetics, 206(4), 2085–2103.
https://doi.org/10.1534/genetics.116.198424

Kalinowski, S. T. (2004). Counting alleles with rarefaction: private
alleles and hierarchical sampling designs. Conservation Genetics, 5(4),
539–543. https://doi.org/10.1023/B:COGE.0000041021.91777.1a

## See also

\[resampling_individuals_CI(), resampling_loci_CI(), calculate_Fst(),
calculate_rarified_Ar()\]
