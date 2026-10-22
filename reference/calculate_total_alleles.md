# Total Allele Counts

This function calculates the number of total and private alleles unique
to each group based on biallelic gentoype data. It computes the count of
private alleles and total alleles for each group and aggregates results
into a data table.

## Usage

``` r
calculate_total_alleles(gt, group_variable, maf)
```

## Arguments

- gt:

  Genotype matrix where individuals are rows and loci are columns, coded
  as 0=aa, 1=aA, 2=AA.

- group_variable:

  Vector of genetic group IDs for each individual (same length as number
  of rows in gt). Can represent sites, populations etc.

- maf:

  Minimum minor allele frequency filter.

## Value

Data table summarizing private alleles and total alleles for each group
and a total count.
