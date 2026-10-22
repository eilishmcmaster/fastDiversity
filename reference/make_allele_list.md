# Create Allele List by Group

This function creates a list of alleles for each group specified by
\`site_variable\` based on genotype data (\`gt\`). It iterates over
unique groups, calculates alleles for each group, and returns a list
where each element corresponds to a group.

## Usage

``` r
make_allele_list(gt, site_variable)
```

## Arguments

- gt:

  Genotype matrix where rows represent individuals and columns represent
  loci (0, 1, 2 format).

- site_variable:

  Vector specifying group IDs (sites) for each individual.

## Value

List where each element contains a vector of alleles for a group.
