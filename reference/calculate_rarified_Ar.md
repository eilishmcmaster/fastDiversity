# Calculate rarified allelic richness

This function calculates rarified allelic richness across genetic groups
(sites) using rarefaction. Rarified allelic richness estimates the
average number of alleles in a sample while accounting for differences
in sample sizes.

## Usage

``` r
calculate_rarified_Ar(gt, site_variable)
```

## Arguments

- gt:

  Genotype matrix where individuals are rows and loci are columns, coded
  as 0=aa, 1=aA, 2=AA.

- site_variable:

  Vector of group IDs (sites) for each individual (same length as number
  of rows in gt).

## Value

A data frame containing rarified allelic richness of each locus per
genetic group (site).

## Details

References: - Kalinowski, S. T. (2005).
https://www.montana.edu/kalinowski/software/documents/2005_HP-Rare_MolecularEcologyNotes.pdf -
Keenan, K., et al. (2013).
https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210x.12067
