# Calculate Allelic Richness

This function calculates allelic richness (Ar) for each locus based on
genotype data. Allelic richness is defined as the average number of
alleles per individual across all loci.

## Usage

``` r
calculate_Ar(gt)
```

## Arguments

- gt:

  Genotype matrix where individuals are rows and loci are columns, coded
  as 0=aa, 1=aA, 2=AA.

## Value

Average allelic richness across loci.

## Details

Reference: - Kalinowski, S. T. (2005).
https://www.montana.edu/kalinowski/software/documents/2005_HP-Rare_MolecularEcologyNotes.pdf
