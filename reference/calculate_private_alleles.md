# Calculate Private Alleles

This function calculates the number of private alleles unique to each
group based on allele lists (\`allele_list\`). It computes the count of
private alleles and total alleles for each group and aggregates results
into a data table.

## Usage

``` r
calculate_private_alleles(allele_list)
```

## Arguments

- allele_list:

  List where each element contains a vector of alleles for a group.

## Value

Data table summarizing private alleles and total alleles for each group
and a total count.
