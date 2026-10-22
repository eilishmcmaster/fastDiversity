# Calculate individual heterozygosity

Calculate individual heterozygosity

## Usage

``` r
individual_Ho(
  gt,
  genetic_group_variable = NULL,
  minimum_loci = 50,
  maf = 0.05,
  max_missingness = 0.3
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

- minimum_loci:

  (optional) Minimum number of loci required to proceed with analysis
  for each genetic group (genetic groups with fewer loci will be
  skipped).

- maf:

  (optional) Minimum minor allele frequency filter.

- max_missingness:

  (optional) Maximum proportion of missing data allowed per locus (loci
  with higher missingness will be filtered out).

## Value

A named vector containing F for each individual in gt.

## References

Frankham, R., Ballou, J. D., & Briscoe, D. A. (2010). Introduction to
Conservation Genetics (2nd ed.). Cambridge University Press.
https://doi.org/10.1017/CBO9780511809002

## Examples

``` r
data(example_data)
# get inbreeding coefficients for individuals by genetic group
Hos <- individual_Ho(example_gt, example_meta$population, maf=0.05)
```
