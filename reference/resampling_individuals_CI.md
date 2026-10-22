# Calculate global bootstrapped He, Ho, and FIS by resampling individuals

Resampling individuals means randomly selecting individuals (rows) from
the dataset with replacement. This keeps the same genetic markers
(loci), but uses different combinations of individuals in each sample.
It tests how your results might change depending on which individuals
you sampled.

## Usage

``` r
resampling_individuals_CI(gt, boots = 100, resample_n = NULL, CI_alpha = 0.05)
```

## Arguments

- gt:

  Genotype matrix where individuals are rows and loci are columns, coded
  as 0=aa, 1=aA, 2=AA.

- boots:

  Number of bootstrap replicates

- resample_n:

  Number of individuals to resample with replacement (defaults to total
  individuals)

- CI_alpha:

  Alpha level for confidence intervals (e.g., 0.05 for 95% CI)

## Value

Named vector of confidence intervals for global FIS, HE, HO

## Details

Bootstrap replicates are drawn with R's random number generator. Call
[`set.seed()`](https://rdrr.io/r/base/Random.html) before this function
for reproducible confidence intervals.
