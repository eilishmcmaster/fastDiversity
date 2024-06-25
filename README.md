fastDiversity: A simple and efficient population genetics analysis package.
=================================================================================

This is the repository for the latest version of the ```fastDiversity``` package. It contains the latest bug fixes, as well as the newest functionality.

This package was created to address shortcomings in existing methods, which often necessitate specific data structures and file types, and are occasionally no longer maintained. It offers essential functions for basic population genomic analysis in a format compatible with any SNP data in biallelic 0,1,2 (aa, aA, AA) format (e.g., DArTseq), thereby saving time compared to alternative approaches.

Installing fastDiversity
=======================
Development version
-------------------

```s
# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install the package
devtools::install("https://github.com/eilishmcmaster/fastDiversity")
```

Example use
=======================
Calculate private and total alleles by specified group (in this case, population):
```
allele_list_population <- make_allele_list(example_gt, example_meta$population)
ggvenn::ggvenn(allele_list_population)
private_total_alleles_population  <- calculate_private_alleles(allele_list_population)
```

Calculate heterozygosity, inbreeding coefficient (Fis), mean allelic richness (Ar), and rarefied allelic richness (rAr):
```
basicstats <- faststats(example_gt, example_meta$population, example_meta$site, minimum_n=3, minimum_loci=50, maf=0.05, max_missingness=0.3)
```


References
=======================

Kalinowski, S. T. (2004). hp‐rare 1.0: a computer program for performing rarefaction on measures of allelic richness. Molecular Ecology Notes, 5(1), 187-189. https://doi.org/10.1111/j.1471-8286.2004.00845.x

Keenan, K. G., McGinnity, P., Cross, T. F., Crozier, W. W., & Prodöhl, P. A. (2013). Diversity: an r package for the estimation and exploration of population genetics parameters and their associated errors. Methods in Ecology and Evolution, 4(8), 782-788. https://doi.org/10.1111/2041-210x.12067
