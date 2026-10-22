# Calculate Expected Heterozygosity per Locus (Hes)

This function calculates expected heterozygosity (Hes) for each locus
based on genotype data. Hes is the expected proportion of heterozygotes
at each locus assuming HW equilibrium.

## Usage

``` r
calculate_Hes(gt)
```

## Arguments

- gt:

  Genotype matrix where individuals are rows and loci are columns, coded
  as 0=aa, 1=aA, 2=AA.

## Value

Vector of expected heterozygosity values (Hes) for each locus.
