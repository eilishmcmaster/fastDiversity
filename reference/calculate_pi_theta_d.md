# Calculate Tajima's D, pi, and Watterson's theta

Calculates Tajima's D, nucleotide diversity (pi), and Watterson's theta
from a biallelic genotype matrix, following the approach of Korunes &
Samuk (2021) to account for missing data.

Monomorphic sites are removed before calculation, and `raw_pi` and
`watterson_theta` are sums over segregating sites. They are not per-site
(per-nucleotide) values, which would require the number of invariant
sites sequenced; this information is unavailable for many
reduced-representation datasets, and monomorphic SNP loci alone do not
provide it. Tajima's D does not depend on invariant sites.

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

  Number of segregating sites

- raw_pi:

  Nucleotide diversity (pi) summed over segregating sites

- watterson_theta:

  Watterson's theta summed over segregating sites

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
