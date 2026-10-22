# Calculate pi, theta, and D statistics

This function calculates pi, Watterson's theta, and Tajima's D for
specified sites within genetic groups following the methods in Pixy,
which accounts for missing data to produce unbiased estimates. For
further details see
https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13326 and
https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.14104 . Data is
filtered based on genetic group (minimum number of loci, minimum minor
allele frequency, and maximum missingness per locus) and minimum number
of individuals per site.

## Usage

``` r
pi_theta_d_stats(
  gt,
  genetic_group_variable,
  site_variable,
  n = NULL,
  minimum_n = 5,
  minimum_loci = 5,
  max_missingness = 0.1,
  maf = 0,
  allele_count_min = 0
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

- n:

  Number of individuals to use for calculations per site. If left as
  NULL, all individuals are used.

- minimum_n:

  Minimum number of individuals required per site (sites with fewer
  individuals will be filtered out).

- minimum_loci:

  Minimum number of loci required to proceed with analysis for each
  genetic group (genetic groups with fewer loci will be skipped).

- max_missingness:

  Maximum proportion of missing data allowed per locus (loci with higher
  missingness will be filtered out).

- maf:

  Minimum minor allele frequency filter.

- allele_count_min:

  Minimum minor allele count to keep a locus, can be used to exclude
  singletons.

## Value

A data frame with the following columns:

- group:

  Genetic group ID (e.g., species, subspecies).

- site:

  Site ID corresponding to the sample location.

- loci:

  Number of loci used in the analysis for that group-site combination.

- n:

  Number of individuals sampled at the site for the given group.

- pi:

  Nucleotide diversity summed over segregating sites.

- S:

  Number of segregating sites.

- theta:

  Watterson's theta summed over segregating sites.

- D:

  Tajima's D statistic.

- real_missingness:

  Proportion of missing data in the filtered genotype matrix.

## Details

**Monomorphic and invariant sites.** Monomorphic sites are excluded, and
`pi` and `theta` are sums over segregating sites rather than per-site
(per-nucleotide) values. Per-site estimates require the number of
invariant sites sequenced, which is unavailable for many
reduced-representation datasets (such as DArTseq SNP reports);
monomorphic SNP loci alone do not provide it. `pi` and `theta` therefore
scale with the number of loci analysed, and are comparable among sites
analysed with the same loci, not with genome-wide nucleotide diversity
estimates. Tajima's D does not depend on invariant sites. See
[`calculate_pi_theta_d`](https://eilishmcmaster.github.io/fastDiversity/reference/calculate_pi_theta_d.md).

## See also

\[resampling_individuals_CI(), resampling_loci_CI()\]
