---
title: 'fastDiversity: an R package for population genetic diversity analysis of multi-species SNP datasets'
tags:
  - R
  - population genetics
  - conservation genetics
  - SNP analysis
  - bioinformatics
authors:
  - name: Eilish S. McMaster
    corresponding: true
    orcid: 0000-0002-7415-8690
    affiliation: 1
affiliations:
  - name: School of Life and Environmental Sciences, University of Sydney, Australia
    index: 1
date: 18 May 2026
bibliography: paper.bib
---

# Summary

Genetic diversity underpins evolutionary potential, adaptive capacity, and population persistence, making its quantification central to conservation genetics and population biology [@frankham2010]. Reduced-representation sequencing approaches including DArTseq [@sansaloni2011] and RADseq [@baird2008] now routinely generate large biallelic single nucleotide polymorphism (SNP) datasets spanning hundreds of individuals and thousands of loci. Increasingly, these datasets include multiple closely related species or divergent lineages analysed within a comparative framework.

fastDiversity is an R package for rapid calculation of population genetic diversity statistics from biallelic SNP dosage matrices. The package operates directly on standard numeric genotype matrices encoded as 0, 1, 2, and NA, avoiding conversion to specialised genetic data classes. fastDiversity implements allele-sharing and private allele analyses; site- and population-level heterozygosity and inbreeding statistics [@nei1978; @hartl2007]; population-specific F~ST~ [@weir2017] and G~ST~ [@nei1973]; rarefied allelic richness [@kalinowski2004]; nucleotide diversity [@nei1979], Watterson’s theta [@watterson1975], and Tajima’s D [@tajima1989], implemented following pixy [@korunes2021]; individual multilocus heterozygosity and inbreeding coefficients [@frankham2010]; and Hardy–Weinberg equilibrium testing [@hardyweinberg], with optional bootstrap confidence intervals.

A defining feature of the package is that locus-level filtering by minor allele frequency, missingness, and allele count is performed independently within user-defined genetic groups. This enables biologically appropriate analysis of multi-species datasets within a single workflow, avoiding artefacts introduced by pooled filtering across reproductively isolated groups.

fastDiversity has been under active development since 2024 and the software is archived on Zenodo (doi:10.5281/zenodo.12741582).

# Statement of need

Several R packages provide population genetic analysis tools for biallelic SNP datasets, including hierfstat [@goudet2005], diveRsity [@keenan2013], adegenet [@jombart2008], dartR [@gruber2018], and poppr [@kamvar2014], collectively covering diversity statistics, ordination, and population structure inference. These tools are widely used but generally apply locus-level filtering globally across the pooled dataset, or require datasets to be manually split by species prior to analysis.

This presents a recurring problem for multi-species datasets. When MAF thresholds, missingness filters, or expected heterozygosity are calculated across pooled samples from divergent taxa, allele frequencies reflect a mixture of reproductively isolated gene pools rather than any biologically meaningful unit. As a result, loci that are polymorphic within a species may be incorrectly removed by pooled MAF filters, while loci that are fixed within species but differentiated between them can inflate diversity or F~IS~ estimates when treated as within-population variation. Differences in locus presence and missingness across taxa further bias locus retention under standard workflows (Table 1). Species-specific filtering avoids these artefacts by evaluating loci within genetic groups, preserving within-taxon signal while enabling comparable summaries across taxa.

Beyond filtering correctness, performing independent analyses per genetic group simultaneously, rather than through repeated manual subsetting, reduces friction and allows comparative diversity analyses to sit alongside other standard steps such as dimensionality reduction, ancestry estimation, and pairwise F~ST~ within a single analysis script.

fastDiversity was developed to address this gap. It applies filtering and diversity estimation independently within user-defined genetic groups while returning results in a unified tidy output, enabling scalable multi-species analyses without repetitive manual re-analysis. The package operates directly on standard SNP dosage matrices rather than requiring conversion to specialised object classes such as genind or genlight, lowering the barrier for users working within existing pipelines.

![**Table 1. Effects of pooled vs species-specific filtering in multi-species SNP datasets.** Example genotype matrix showing how locus informativeness differs across three species (A–C) under pooled filtering. Loci may be informative across all species (L1), polymorphic only within one species and potentially removed by pooled MAF filters (L2), alternatively fixed within species inflating pooled summary statistics such as FIS (L3), informative in one species but fixed in others  and potentially removed by pooled MAF filters (L4), or present only in a subset of species and excluded under pooled missingness filtering (L5). Species-specific filtering retains different informative loci per taxon (A: L1, L4; B: L1, L5; C: L1, L2), illustrating how pooled filtering can distort comparative SNP datasets.](table1.png){width=70%}



# Software design and functionality

fastDiversity is designed around two principles: biologically appropriate filtering and minimal data-handling overhead.

The package accepts standard biallelic SNP dosage matrices produced by common variant-calling pipelines, together with plain metadata vectors describing sample groupings. This design allows integration into existing R analysis pipelines without object conversion or restructuring.

All locus-level quality filters, including missingness and minor allele frequency thresholds, are applied independently within each user-defined genetic group. This filtering framework propagates throughout the package and ensures that downstream diversity statistics are calculated only from loci informative within the relevant biological unit.

Because most SNP-calling pipelines report only polymorphic sites, the documentation distinguishes statistics by their dependence on monomorphic sites. Ratio statistics (F~IS~, G~ST~, F~ST~) and Tajima's D are unaffected by them. Heterozygosity and allelic richness are averages over the loci analysed, so they are comparable among sites within a genetic group but are not genome-wide estimates. Nucleotide diversity and Watterson's θ are calculated from segregating sites only and reported as sums rather than per-site values, because the number of invariant sites sequenced is unavailable for many reduced-representation datasets.

Implemented analyses span multiple biological scales, including allele-level summaries such as allele sharing and private alleles; population-level diversity statistics including observed and expected heterozygosity, inbreeding coefficients, rarefied allelic richness [@kalinowski2004], and differentiation as population-specific F~ST~ [@weir2017], the per-population counterpart of Weir and Cockerham's θ [@weir1984], and G~ST~ [@nei1973]; sequence diversity statistics including nucleotide diversity, Watterson’s theta, and Tajima’s D using missing-data-aware approaches [@korunes2021];
individual-level summaries including multilocus heterozygosity and Wright’s inbreeding coefficient. Where appropriate, bootstrap and resampling procedures adapted from @goudet2005 and @keenan2013 are implemented to estimate confidence intervals and standardise comparisons among uneven sample sizes.

Summary functions return tidy data frames [@wickham2014], with one row per genetic group and site and one column per statistic, which can be passed directly to standard R plotting and data-manipulation tools.

# State of the field

hierfstat [@goudet2005], diveRsity [@keenan2013], adegenet [@jombart2008], poppr [@kamvar2014], dartR [@gruber2018], and tidypopgen [@carter2025] provide established implementations of diversity and differentiation statistics, but each retains a single set of loci for all individuals in a dataset. dartR can take populations into account when filtering by minor allele frequency, but the retained loci are still shared across all individuals, and locus selection in tidypopgen likewise applies to the whole dataset. Applying the same filters independently within each species therefore requires splitting the data and repeating the analysis for each group. fastDiversity applies the filters within every genetic group in a single call and returns the results in one table, complementing these packages for comparative multi-species SNP datasets.

# Research impact statement

fastDiversity has been used in more than 30 confidential threatened-species assessments for the New South Wales Government Saving our Species program [@doyle2025; @cascini2025], where comparable diversity estimates across related taxa inform conservation management. It has also supported peer-reviewed conservation genomic studies of threatened Australian plants, including multi-species comparisons of diversity and population structure in *Homoranthus* (Myrtaceae) [@mcmaster2025homoranthus], subspecies delimitation in the critically endangered *Zieria obcordata* (Rutaceae) [@mcmaster2025zieria], and evaluation of genetic rescue risks for Australia’s rarest *Eucalyptus* [@mcmaster2026eucalyptus], as well as a comparison of kinship estimation methods for reduced-representation SNP data [@mcmaster2025kinship].

# Acknowledgements

Development of fastDiversity benefited from discussions with collaborators in conservation genomics and population genetics at the University of Sydney, the Botanic Gardens of Sydney, and Queensland Herbarium.

# AI usage disclosure

Claude Sonnet 4.6 (Anthropic) and ChatGPT (OpenAI GPT-5.5) were used to assist with editing and restructuring documentation and manuscript text based on author-provided drafts and outlines. All scientific content, software design decisions, implementations, analyses, and interpretations were developed, verified, and approved by the author.

# References
