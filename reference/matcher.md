# Identify Alleles Matching Genotype

This function identifies alleles at specified loci that match the
genotype data (\`gt\`). It checks the genotype matrix for each locus and
retrieves corresponding alleles from a provided loci data frame
(\`loci\`).

## Usage

``` r
matcher(gt, loci)
```

## Arguments

- gt:

  Genotype matrix where the first row represents loci and subsequent
  rows are genotype data coded as 0=aa, 1=aA, 2=AA.

- loci:

  Data frame containing locus information with columns 'loci',
  'allele1', and 'allele2'.

## Value

Vector of alleles matching the genotype at each locus in \`gt\`.
