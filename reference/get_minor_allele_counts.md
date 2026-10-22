# Get Minor Allele Counts

This function calculates the minor allele counts for each locus based on
genotype data. The minor allele count is defined as the sum of the less
common allele at each locus.

## Usage

``` r
get_minor_allele_counts(gt)
```

## Arguments

- gt:

  Genotype matrix where individuals are rows and loci are columns, coded
  as 0=aa, 1=aA, 2=AA.

## Value

Vector of minor allele counts for each locus.
