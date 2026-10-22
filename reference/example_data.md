# Example SNP genotype data and sample metadata

DArTseq SNP genotypes for 257 individuals of an unpublished plant
species collected in eastern Australia, with matching sample metadata.
The species identity is withheld pending publication.
`data(example_data)` loads two objects, `example_gt` and `example_meta`.

## Usage

``` r
data(example_data)
```

## Format

`example_gt`: a data frame with 257 rows (individuals) and 1906 columns
(SNP loci), with genotypes coded as 0 (homozygous reference, aa), 1
(heterozygous, aA), 2 (homozygous alternate, AA) and `NA` (missing;
about 13% of genotypes). Row names are sample IDs and column names are
DArTseq locus IDs. All loci are polymorphic across the full dataset;
monomorphic sites are not included.

`example_meta`: a data frame with 257 rows, in the same order as the
rows of `example_gt`, and 3 columns:

- sample:

  Sample ID, matching `rownames(example_gt)`.

- site:

  Collection site (55 sites).

- population:

  Genetic group assignment (integer, 5 groups).

## Source

Unpublished DArTseq data provided by the package author.

## Examples

``` r
data(example_data)
dim(example_gt)
#> [1]  257 1906
example_gt[1:5, 1:5]
#>            X26095427.F.0.51.G.T.51.G.T X26097233.F.0.11.C.A.11.C.A
#> NSW1008199                           0                           2
#> NSW1008169                           0                           2
#> NSW1008189                           0                           2
#> NSW1008229                           0                           2
#> NSW1008224                           0                           2
#>            X26097154.F.0.27.G.C.27.G.C X26099383.F.0.35.A.G.35.A.G
#> NSW1008199                           0                           2
#> NSW1008169                           0                           2
#> NSW1008189                           0                           2
#> NSW1008229                           0                           2
#> NSW1008224                           0                           2
#>            X26101203.F.0.30.A.C.30.A.C
#> NSW1008199                           0
#> NSW1008169                           0
#> NSW1008189                           0
#> NSW1008229                           0
#> NSW1008224                           0
head(example_meta)
#>       sample         site population
#> 1 NSW1008199 Wambam_Creek          5
#> 2 NSW1008169 Wambam_Creek          5
#> 3 NSW1008189 Wambam_Creek          5
#> 4 NSW1008229 Wambam_Creek          5
#> 5 NSW1008224 Wambam_Creek          5
#> 6 NSW1009574 Wambam_Creek          5
table(example_meta$population)
#> 
#>   1   2   3   4   5 
#> 100  23  46  31  57 
```
