# Calculate individual inbreeding coefficient (Wright's F)

This function calculates Wright's F (individual inbreeding coefficient),
the probability that two alleles at a locus in an individual are
identical by descent. Calculation from Frankham, Ballou, & Briscoe
(2010, p 266) F= 1-(Ho/He).

## Usage

``` r
individual_F(
  gt,
  genetic_group_variable = NULL,
  site_variable = NULL,
  minimum_loci = 50,
  maf = 0.05,
  max_missingness = 0.3,
  minimum_n = 3,
  table_out = FALSE
)
```

## Arguments

- gt:

  Genotype matrix where individuals are rows and loci are columns, coded
  as 0=aa, 1=aA, 2=AA.

- genetic_group_variable:

  (optional) Vector of genetic group IDs for each individual (same
  length as number of rows in gt). Can represent species, subspecies, or
  other biologically relevant structure groups. If not included, the
  whole gt will be treated as one group.

- site_variable:

  (optional) Vector of site IDs for each individual. If supplied,
  expected heterozygosity will be calculated within sites, which will
  probably result in F being lower due to the Wahlund effect.

- minimum_loci:

  (optional) Minimum number of loci required to proceed with analysis
  for each genetic group (genetic groups with fewer loci will be
  skipped).

- maf:

  (optional) Minimum minor allele frequency filter.

- max_missingness:

  (optional) Maximum proportion of missing data allowed per locus (loci
  with higher missingness will be filtered out).

- minimum_n:

  (optional) Minimum number of individuals required in a group or site
  for F to be calculated (individuals in smaller groups or sites return
  NA).

- table_out:

  (optional) If TRUE, return a data frame with columns individual, F,
  group and site instead of a named vector.

## Value

A named vector containing F for each individual in gt, or a data frame
if `table_out = TRUE`.

## References

Frankham, R., Ballou, J. D., & Briscoe, D. A. (2010). Introduction to
Conservation Genetics (2nd ed.). Cambridge University Press.
https://doi.org/10.1017/CBO9780511809002

## Examples

``` r
data(example_data)
# get inbreeding coefficients for individuals by genetic group
Fs <- individual_F(example_gt, example_meta$population)
```
