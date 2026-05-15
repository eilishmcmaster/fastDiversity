[![DOI](https://zenodo.org/badge/797486652.svg)](https://zenodo.org/doi/10.5281/zenodo.12741582)

# fastDiversity

**fastDiversity** is an R package for fast, flexible calculation of population genetic diversity statistics from biallelic SNP data.

Many existing tools require specific file formats or data structures, and some are no longer actively maintained. fastDiversity is designed to work directly with any SNP dataset encoded in standard biallelic dosage format (`0` (homozygous reference, *aa*), `1` (heterozygous, *aA*), `2` (homozygous alternate, *AA*), with `NA` for missing data) including the direct output of DArTseq and most SNP-calling pipelines.

A key design goal is support for datasets spanning **multiple species or genetic groups**, where locus filtering and statistics need to be calculated independently per group. This is a common scenario in conservation genomics and comparative population studies that is not well served by existing tools.

## Functions

| Function | Description |
|----------|-------------|
| `make_allele_list()` | Returns a named list of allele sets per group, for allele sharing analysis |
| `calculate_private_alleles()` | Counts private (group-exclusive) and total alleles per group |
| `faststats()` | Calculates Ho, He, uHe, Fis, uFis, Fst, and Ar per site, with per-group locus filtering, optional bootstrapped confidence intervals, and optional Hardy-Weinberg equilibrium testing |

## Installation

```r
# Install remotes if needed
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

remotes::install_github("eilishmcmaster/fastDiversity")
```

## Quick start

```r
library(fastDiversity)

data(example_gt)    # genotype matrix: samples × loci, values 0/1/2/NA
data(example_meta)  # sample metadata with population and site columns
```

**Allele sharing and private alleles:**

```r
allele_list <- make_allele_list(example_gt, example_meta$population)

# Visualise overlap between populations
ggvenn::ggvenn(allele_list)

# Count private and total alleles per population
calculate_private_alleles(allele_list)
```

**Diversity statistics:**

```r
basicstats <- faststats(
  gt                     = example_gt,
  genetic_group_variable = example_meta$population,
  site_variable          = example_meta$site,
  minimum_n              = 3,
  minimum_loci           = 50,
  maf                    = 0.05,
  max_missingness        = 0.3
)
```

Returns one row per genetic group × site combination, with columns for Ho, He, uHe, Fis, uFis, Fst, Ar, and the number of individuals and loci used.

**With bootstrapped confidence intervals:**

```r
basicstats_ci <- faststats(
  gt                     = example_gt,
  genetic_group_variable = example_meta$population,
  site_variable          = example_meta$site,
  minimum_n              = 3,
  minimum_loci           = 50,
  maf                    = 0.05,
  max_missingness        = 0.3,
  get_CI_loci_resampling = TRUE,   # resample loci
  boots                  = 100,
  CI_alpha               = 0.05
)
```
**Individual-level statistics:**

```r
# Individual multilocus heterozygosity (MLH)
ind_MLH_ho <- individual_Ho(
  example_gt,
  genetic_group_variable = example_meta$population,
  maf                    = 0.05,
  max_missingness        = 0.3
)

# Individual inbreeding coefficient (Wright's F = 1 - Ho/He)
# He is calculated at the genetic group level, avoiding Wahlund effect inflation
ind_F <- individual_F(
  example_gt,
  genetic_group_variable = example_meta$population,
  maf                    = 0.05,
  max_missingness        = 0.3
)

plot(ind_MLH_ho, ind_F)
```

Both functions apply MAF and missingness filtering per genetic group before
calculation, so individual estimates are based on the same loci used in
population-level analyses.


For a full worked example with plots, see `vignette("fastDiversity-intro")`.

## References

Keenan, K. G., McGinnity, P., Cross, T. F., Crozier, W. W., & Prodöhl, P. A. (2013). diveRsity: An R package for the estimation and exploration of population genetics parameters and their associated errors. *Methods in Ecology and Evolution*, 4(8), 782–788. https://doi.org/10.1111/2041-210x.12067

Nei, M. (1978). Estimation of average heterozygosity and genetic distance from a small number of individuals. *Genetics*, 89(3), 583–590. https://doi.org/10.1093/genetics/89.3.583

Weir, B. S., & Cockerham, C. C. (1984). Estimating F-statistics for the analysis of population structure. *Evolution*, 38(6), 1358–1370. https://doi.org/10.1111/j.1558-5646.1984.tb05657.x