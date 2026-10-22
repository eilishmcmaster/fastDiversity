# Test for Hardy-Weinberg Equilibrium (HWE) across multiple loci using HardyWeinberg

This function performs Hardy-Weinberg equilibrium tests on multiple
genetic loci using the HWExact from the HardyWeinberg package. It
returns the Fisher combined p-value of the LLR test, the mean U-score,
and individual test results for each locus.

## Usage

``` r
MC_HWEtest_of_loci(gt)
```

## Arguments

- gt:

  A data frame or matrix where each column corresponds to genotypes at a
  single locus. Genotypes must be coded as 0 (AA), 1 (AB), or 2 (BB).

## Value

A list containing:

- HW_fisher_p:

  Fisher combined p-value from the HWExact test across loci.

- loci_HWX_results:

  A data table with individual p-values and locus names.
