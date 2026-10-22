# Get Minor Allele Frequencies

This function calculates the minor allele frequencies for each locus
based on genotype data. The minor allele frequency is defined as the
frequency of the less common allele at each locus.

## Usage

``` r
get_minor_allele_frequencies(gt)
```

## Arguments

- gt:

  Genotype matrix where individuals are rows and loci are columns, coded
  as 0=aa, 1=aA, 2=AA.

## Value

Vector of minor allele frequencies for each locus.
