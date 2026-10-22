# Package index

## Diversity statistics by group and site

Main wrappers that filter loci per genetic group and summarise
statistics per site.

- [`faststats()`](https://eilishmcmaster.github.io/fastDiversity/reference/faststats.md)
  : Calculate basic diversity statistics
- [`pi_theta_d_stats()`](https://eilishmcmaster.github.io/fastDiversity/reference/pi_theta_d_stats.md)
  : Calculate pi, theta, and D statistics

## Individual-level statistics

- [`individual_F()`](https://eilishmcmaster.github.io/fastDiversity/reference/individual_F.md)
  : Calculate individual inbreeding coefficient (Wright's F)
- [`individual_Ho()`](https://eilishmcmaster.github.io/fastDiversity/reference/individual_Ho.md)
  : Calculate individual heterozygosity

## Allele sharing and private alleles

- [`make_allele_list()`](https://eilishmcmaster.github.io/fastDiversity/reference/make_allele_list.md)
  : Create Allele List by Group
- [`calculate_private_alleles()`](https://eilishmcmaster.github.io/fastDiversity/reference/calculate_private_alleles.md)
  : Calculate Private Alleles
- [`calculate_total_alleles()`](https://eilishmcmaster.github.io/fastDiversity/reference/calculate_total_alleles.md)
  : Total Allele Counts
- [`matcher()`](https://eilishmcmaster.github.io/fastDiversity/reference/matcher.md)
  : Identify Alleles Matching Genotype

## Confidence intervals

- [`resampling_loci_CI()`](https://eilishmcmaster.github.io/fastDiversity/reference/resampling_loci_CI.md)
  : Calculate global bootstrapped He, Ho, and FIS by resampling loci (as
  in hierfstat::boot.ppfis or arlequin)
- [`resampling_individuals_CI()`](https://eilishmcmaster.github.io/fastDiversity/reference/resampling_individuals_CI.md)
  : Calculate global bootstrapped He, Ho, and FIS by resampling
  individuals

## Hardy-Weinberg equilibrium

- [`HWEtest_of_loci()`](https://eilishmcmaster.github.io/fastDiversity/reference/HWEtest_of_loci.md)
  : Test for Hardy-Weinberg Equilibrium (HWE) across multiple loci using
  HWxtest
- [`MC_HWEtest_of_loci()`](https://eilishmcmaster.github.io/fastDiversity/reference/MC_HWEtest_of_loci.md)
  : Test for Hardy-Weinberg Equilibrium (HWE) across multiple loci using
  HardyWeinberg

## Building blocks

Lower-level functions used by the wrappers above.

- [`calculate_Ho()`](https://eilishmcmaster.github.io/fastDiversity/reference/calculate_Ho.md)
  : Calculate Observed Heterozygosity (Ho)
- [`calculate_Hes()`](https://eilishmcmaster.github.io/fastDiversity/reference/calculate_Hes.md)
  : Calculate Expected Heterozygosity per Locus (Hes)
- [`calculate_uHe()`](https://eilishmcmaster.github.io/fastDiversity/reference/calculate_uHe.md)
  : Calculate Unbiased Expected Heterozygosity (uHe)
- [`calculate_Fst()`](https://eilishmcmaster.github.io/fastDiversity/reference/calculate_Fst.md)
  : Calculate population-specific Fst (Weir & Goudet 2017)
- [`calculate_Ar()`](https://eilishmcmaster.github.io/fastDiversity/reference/calculate_Ar.md)
  : Calculate Allelic Richness (unrarefied)
- [`calculate_rarified_Ar()`](https://eilishmcmaster.github.io/fastDiversity/reference/calculate_rarified_Ar.md)
  : Calculate rarified allelic richness
- [`calculate_pi_theta_d()`](https://eilishmcmaster.github.io/fastDiversity/reference/calculate_pi_theta_d.md)
  : Calculate Tajima's D, pi, and Watterson's theta
- [`get_minor_allele_counts()`](https://eilishmcmaster.github.io/fastDiversity/reference/get_minor_allele_counts.md)
  : Get Minor Allele Counts
- [`get_minor_allele_frequencies()`](https://eilishmcmaster.github.io/fastDiversity/reference/get_minor_allele_frequencies.md)
  : Get Minor Allele Frequencies

## Data

- [`example_data`](https://eilishmcmaster.github.io/fastDiversity/reference/example_data.md)
  : Example SNP genotype data and sample metadata

## Package

- [`fastDiversity`](https://eilishmcmaster.github.io/fastDiversity/reference/fastDiversity-package.md)
  [`fastDiversity-package`](https://eilishmcmaster.github.io/fastDiversity/reference/fastDiversity-package.md)
  : fastDiversity: Population Genetics Analysis
