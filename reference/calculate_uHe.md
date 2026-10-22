# Calculate Unbiased Expected Heterozygosity (uHe)

This function calculates unbiased expected heterozygosity (uHe) based on
the number of individuals (ns) and observed heterozygosity (Hes) for
each genetic group or locus. Unbiased expected heterozygosity corrects
for sample size biases and provides an estimate of heterozygosity
accounting for the number of alleles and individuals.

## Usage

``` r
calculate_uHe(ns, Hes)
```

## Arguments

- ns:

  Vector containing the number of individuals sampled for each genetic
  group or locus.

- Hes:

  Vector containing observed heterozygosity values for each genetic
  group or locus.

## Value

Average unbiased expected heterozygosity (uHe) across genetic groups or
loci.
