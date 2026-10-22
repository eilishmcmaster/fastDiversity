# Calculate Allelic Richness (unrarefied)

This function calculates allelic richness (Ar) as the number of distinct
alleles observed at each locus (1 or 2 for biallelic SNPs), averaged
across loci. This value is not corrected for sample size, so it should
only be compared between samples of similar size. Use
[`calculate_rarified_Ar`](https://eilishmcmaster.github.io/fastDiversity/reference/calculate_rarified_Ar.md)
for allelic richness standardised to a common sample size, as reported
by
[`faststats`](https://eilishmcmaster.github.io/fastDiversity/reference/faststats.md).

## Usage

``` r
calculate_Ar(gt)
```

## Arguments

- gt:

  Genotype matrix where individuals are rows and loci are columns, coded
  as 0=aa, 1=aA, 2=AA.

## Value

Mean number of alleles observed per locus. Loci with no genotypes are
ignored.

## See also

[`calculate_rarified_Ar`](https://eilishmcmaster.github.io/fastDiversity/reference/calculate_rarified_Ar.md)
