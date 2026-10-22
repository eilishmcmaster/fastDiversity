# Test for Hardy-Weinberg Equilibrium (HWE) across multiple loci using HWxtest

**This function is defunct.** This function performs Hardy-Weinberg
equilibrium tests on multiple genetic loci using the hwx.test from the
HWxtest package. It returns the Fisher combined p-value of the LLR test,
the mean U-score, and individual test results for each locus.

## Usage

``` r
HWEtest_of_loci(gt)
```

## Arguments

- gt:

  A data frame or matrix where each column corresponds to genotypes at a
  single locus. Genotypes must be coded as 0 (AA), 1 (AB), or 2 (BB).

## Value

A named vector containing:
