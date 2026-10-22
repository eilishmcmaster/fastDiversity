<!-- badges: start -->
[![DOI](https://zenodo.org/badge/797486652.svg)](https://zenodo.org/doi/10.5281/zenodo.12741582)
[![R-CMD-check](https://github.com/eilishmcmaster/fastDiversity/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/eilishmcmaster/fastDiversity/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

# fastDiversity

**Documentation:** <https://eilishmcmaster.github.io/fastDiversity/>

**fastDiversity** is an R package for fast, flexible calculation of population genetic diversity statistics from biallelic SNP data.

Many existing tools require specific file formats or data structures, and some are no longer actively maintained. fastDiversity is designed to work directly with any SNP dataset encoded in standard biallelic dosage format (`0` (homozygous reference, *aa*), `1` (heterozygous, *aA*), `2` (homozygous alternate, *AA*), with `NA` for missing data), the direct output of most SNP-calling pipelines.

A key design goal is support for datasets spanning **multiple species or genetic groups**, where locus filtering and statistics need to be calculated independently per group. This is a common scenario in conservation genomics and comparative population studies that is not well served by existing tools.

## Key functions

| Function | Description |
|----------|-------------|
| `make_allele_list()` | Returns a named list of allele sets per group, for allele sharing analysis |
| `calculate_private_alleles()` | Counts private (group-exclusive) and total alleles per group |
| `faststats()` | Calculates Ho, He, uHe, Fis, uFis, Gst, Fst, and rarefied Ar per site, with per-group locus filtering, optional bootstrapped confidence intervals, and optional Hardy-Weinberg equilibrium testing |
| `pi_theta_d_stats()` | Calculates pi, theta, S, and D with per-group locus filtering |

## Installation

```r
# Install remotes if needed
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

remotes::install_github("eilishmcmaster/fastDiversity")

library(fastDiversity)
```

The tutorial vignette is available online at
<https://eilishmcmaster.github.io/fastDiversity/articles/fastDiversity-intro.html>.

To build the vignette locally, install its additional dependencies first:

```r
install.packages(c("knitr", "rmarkdown", "ggplot2", "ggthemes", "tidyr",
                   "ggvenn", "ggVennDiagram"))
remotes::install_github("eilishmcmaster/fastDiversity", build_vignettes = TRUE)
browseVignettes("fastDiversity")
```

## Quick start

```r
library(fastDiversity)

data(example_data)  # loads two objects:
                    #   example_gt:   genotype matrix, samples × loci, values 0/1/2/NA
                    #   example_meta: sample metadata with population and site columns
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

Returns one row per genetic group × site combination, with columns for Ho, He, uHe, Fis, uFis, Gst (Nei 1973), population-specific Fst (Weir & Goudet 2017), rarefied allelic richness (Ar), and the number of individuals and loci used.

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

***Neutrality statistics:***

```r
# Use all available individuals per population
D_n_all <- pi_theta_d_stats(
  example_gt,
  genetic_group_variable = example_meta$population,
  site_variable          = example_meta$population,
  n                      = NULL, # if set to NULL, all individuals are used, but downsampling can be specified
  minimum_n              = 5, 
  minimum_loci           = 5,
  max_missingness        = 0.1,
  maf                    = 0,
  allele_count_min       = 0
)
```
Returns a table containing nucleotide diversity (pi), segregating sites (S), 
Wu & Watterson's theta (W), and Tajima's D (D) per site_variable grouping. 


For a full worked example with plots, see the
[tutorial vignette](https://eilishmcmaster.github.io/fastDiversity/articles/fastDiversity-intro.html).

## Development and testing

The package is tested with [testthat](https://testthat.r-lib.org/). From a
clone of the repository, run:

```r
devtools::test()   # run the test suite
devtools::check()  # full R CMD check, including the tests
```

R CMD check runs automatically on macOS, Windows, and Linux via GitHub Actions
for every push and pull request to `main`.

## Contributing and support

- **Bug reports and questions:** please open an
  [issue](https://github.com/eilishmcmaster/fastDiversity/issues), including a
  minimal reproducible example and the output of `sessionInfo()`.
- **Contributions:** pull requests are welcome. See
  [CONTRIBUTING.md](.github/CONTRIBUTING.md) for guidelines.

## References

Goudet, J. (2005). Hierfstat, a package for r to compute and test hierarchical F-statistics. 
*Molecular Ecology Notes*, 5(1), 184–186. https://doi.org/10.1111/j.1471-8286.2004.00828.x

Kalinowski, S. T. (2004). hp-rare 1.0: a computer program for performing
rarefaction on measures of allelic richness. *Molecular Ecology Notes*, 5(1),
187–189. https://doi.org/10.1111/j.1471-8286.2004.00845.x

Keenan, K. G., McGinnity, P., Cross, T. F., Crozier, W. W., & Prodöhl, P. A.
(2013). diveRsity: An R package for the estimation and exploration of
population genetics parameters and their associated errors. *Methods in
Ecology and Evolution*, 4(8), 782–788.
https://doi.org/10.1111/2041-210x.12067

Korunes, K. L., & Samuk, K. (2021). pixy: Unbiased estimation of nucleotide diversity 
and divergence in the presence of missing data. *Molecular Ecology Resources*,
21(4), 1359–1368. https://doi.org/10.1111/1755-0998.13326

Nei, M. (1973). Analysis of gene diversity in subdivided populations.
*Proceedings of the National Academy of Sciences*, 70(12), 3321–3323.
https://doi.org/10.1073/pnas.70.12.3321

Nei, M. (1978). Estimation of average heterozygosity and genetic distance from
a small number of individuals. *Genetics*, 89(3), 583–590.
https://doi.org/10.1093/genetics/89.3.583

Weir, B. S., & Goudet, J. (2017). A unified characterization of population
structure and relatedness. *Genetics*, 206(4), 2085–2103.
https://doi.org/10.1534/genetics.116.198424
