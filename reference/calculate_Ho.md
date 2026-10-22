# Calculate Observed Heterozygosity (Ho)

This function calculates observed heterozygosity (Ho) based on genotype
data. Heterozygosity (Ho) is defined as the proportion of heterozygous
individuals among all individuals with genotype data.

## Usage

``` r
calculate_Ho(gt)
```

## Arguments

- gt:

  Genotype matrix where individuals are rows and loci are columns, coded
  as 0=aa, 1=aA, 2=AA.

## Value

Observed heterozygosity (Ho) across all loci.
