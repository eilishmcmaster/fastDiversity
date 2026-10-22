# Calculate population-specific Fst (Weir & Goudet 2017)

This function calculates population-specific Fst for each site using the
allele-matching estimator of Weir & Goudet (2017). For each site, Fst
compares the probability that two alleles drawn from different
individuals within that site are identical (within-site matching) with
the average probability that two alleles drawn from different sites are
identical (between-site matching): \\F\_{ST,i} = (M_i - M_B) / (1 -
M_B)\\.

## Usage

``` r
calculate_Fst(gt, site_variable)
```

## Arguments

- gt:

  Genotype matrix where individuals are rows and loci are columns, coded
  as 0=aa, 1=aA, 2=AA.

- site_variable:

  Vector of site IDs for each individual (same length as number of rows
  in gt). At least two sites are required.

## Value

Named vector of Fst values, one per site.

## Details

This is the population-specific counterpart of Weir & Cockerham's (1984)
\\\theta\\, and gives the same values as the `Fs["Fst", ]` row of
[`hierfstat::fs.dosage()`](https://rdrr.io/pkg/hierfstat/man/fs.dosage.html).
Missing genotypes are handled as in hierfstat, by averaging matching
over the loci genotyped in both individuals of a pair. Values are
relative to the set of sites supplied, so they depend on which sites are
included.

## References

Weir, B. S., & Goudet, J. (2017). A unified characterization of
population structure and relatedness. Genetics, 206(4), 2085–2103.
https://doi.org/10.1534/genetics.116.198424

Weir, B. S., & Cockerham, C. C. (1984). Estimating F-statistics for the
analysis of population structure. Evolution, 38(6), 1358–1370.
https://doi.org/10.2307/2408641

## Examples

``` r
data(example_data)
keep <- example_meta$population == 1
calculate_Fst(example_gt[keep, ], example_meta$site[keep])
#>             Killen_Falls          Big_Scrub Loop1          Big_Scrub Loop2 
#>                      NaN              0.060757857              0.012385405 
#>              Old_Mill_Rd             Rummery_Road       Tungun Rd Rocky Ck 
#>              0.118468010              0.067959552              0.006383656 
#>            Terania_Creek Gibbergunyah Nightcap Tk            Koonyum_Range 
#>              0.080483435              0.018189647              0.173950125 
#>                 Mt_Nardi             Mt_Jerusalem                Blackbutt 
#>              0.067128383              0.167004634              0.036830140 
#>               Airdrop_Rd            Mt_Nothofagus      Best_of_All_Lookout 
#>             -0.003439167              0.087569096              0.067493223 
#>             Lamington_NP              Binna_Burra          Mary_Cairncross 
#>              0.018097713              0.039678375              0.225912684 
```
