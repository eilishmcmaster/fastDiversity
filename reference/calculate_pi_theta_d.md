# Calculate Tajima's D, pi, and Watterson's theta

Calculates Tajima's D, nucleotide diversity (pi), and Watterson's theta
from a biallelic genotype matrix, following the approach of Korunes &
Samuk (2021) to account for missing data. Non-segregating and
monomorphic sites are removed prior to calculation.

## Usage

``` r
calculate_pi_theta_d(gt)
```

## Arguments

- gt:

  Genotype matrix where individuals are rows and loci are columns, coded
  as 0=aa, 1=aA, 2=AA.

## Value

A list with elements:

- tajima_d:

  Tajima's D statistic

- num_sites:

  Number of sites used after filtering

- raw_pi:

  Nucleotide diversity (pi)

- watterson_theta:

  Watterson's theta

- d_stdev:

  Standard deviation of the denominator

## References

Korunes, K. L., & Samuk, K. (2021). pixy: Unbiased estimation of
nucleotide diversity and divergence in the presence of missing data.
Molecular Ecology Resources, 21(4), 1359–1368.
https://doi.org/10.1111/1755-0998.13326

Tajima, F. (1989). Statistical method for testing the neutral mutation
hypothesis by DNA polymorphism. Genetics, 123(3), 585–595.
https://doi.org/10.1093/genetics/123.3.585
