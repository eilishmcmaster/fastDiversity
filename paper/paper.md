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

fastDiversity is an R package for rapid calculation of population genetic diversity statistics from biallelic SNP dosage matrices. The package operates directly on standard numeric genotype matrices encoded as 0, 1, 2, and NA, avoiding conversion to specialised genetic data classes. fastDiversity implements allele-sharing and private allele analyses; site- and population-level heterozygosity and inbreeding statistics; Weir–Cockerham F-statistics [@weir1984]; nucleotide diversity, Watterson’s theta, and Tajima’s D; multilocus heterozygosity; and Hardy–Weinberg equilibrium testing with optional bootstrap confidence intervals.

A defining feature of the package is that locus-level filtering by minor allele frequency, missingness, and allele count is performed independently within user-defined genetic groups. This enables biologically appropriate analysis of multi-species datasets within a single workflow, avoiding artefacts introduced by pooled filtering across reproductively isolated groups.

fastDiversity has been under active development since 2024 and has been applied in conservation genomics studies of threatened Australian plant species and comparative multi-species datasets [@mcmaster2025homoranthus; @mcmaster2025zieria; @mcmaster2025kinship].The software is archived on Zenodo (doi:10.5281/zenodo.12741582).

# Statement of need

Several R packages provide population genetic analysis tools for SNP datasets, including hierfstat [@goudet2005], diveRsity [@keenan2013], adegenet [@jombart2008], dartR [@gruber2018], and poppr [@kamvar2014]. These packages are widely used and provide extensive functionality for diversity analysis, ordination, and population structure inference.

However, multi-species SNP datasets present a recurring practical problem. Existing workflows generally apply locus filtering globally across the pooled dataset or require manual splitting of datasets by species prior to analysis. For comparative datasets spanning multiple species or divergent lineages, pooled filtering can produce biologically misleading estimates. For example, a locus that is monomorphic in one species but polymorphic in another may pass or fail filtering thresholds based on pooled allele frequencies that are representative of neither group individually. Similarly, expected heterozygosity calculated across reproductively isolated taxa lacks clear biological interpretation.

fastDiversity was developed to address this gap by applying filtering and diversity estimation independently within user-defined genetic groups while returning results in a unified tidy output structure. This enables scalable comparative analyses across many species or populations without repetitive manual subsetting and re-analysis.

The package additionally lowers barriers to analysis by operating directly on standard SNP dosage matrices rather than requiring conversion into specialised object classes such as genind or genlight.

# Software design and functionality

fastDiversity is designed around two principles: biologically appropriate filtering and minimal data-handling overhead.

The package accepts standard biallelic SNP dosage matrices produced by common variant-calling pipelines, together with plain metadata vectors describing sample groupings. This design allows integration into existing R analysis pipelines without object conversion or restructuring.

All locus-level quality filters, including missingness and minor allele frequency thresholds, are applied independently within each user-defined genetic group. This filtering framework propagates throughout the package and ensures that downstream diversity statistics are calculated only from loci informative within the relevant biological unit.

Implemented analyses span multiple biological scales, including allele-level summaries such as allele sharing and private alleles; population-level diversity statistics including observed and expected heterozygosity, inbreeding coefficients, and Weir–Cockerham F-statistics [@weir1984]; sequence diversity statistics including nucleotide diversity, Watterson’s theta, and Tajima’s D using missing-data-aware approaches [@korunes2021];
individual-level summaries including multilocus heterozygosity and Wright’s inbreeding coefficient. Where appropriate, bootstrap and resampling procedures adapted from @goudet2005 and @keenan2013 are implemented to estimate confidence intervals and standardise comparisons among uneven sample sizes.

All functions return tidy tabular outputs compatible with standard R visualisation and downstream analysis workflows.

# State of the field

fastDiversity complements rather than replaces existing population genetic software ecosystems. hierfstat [@goudet2005] and diveRsity [@keenan2013] provide established implementations of diversity and differentiation statistics but apply filtering globally and require specialised input structures. adegenet [@jombart2008] and dartR [@gruber2018] provide extensive infrastructure for genomic analyses and DArTseq workflows but are centred on custom object classes and single-species analytical assumptions. poppr [@kamvar2014] supports clonal and mixed-ploidy systems but does not specifically address per-group filtering in comparative multi-species datasets.

fastDiversity occupies a complementary niche focused on lightweight analysis of comparative SNP datasets in which biologically scoped filtering is required across many taxa or lineages simultaneously.

# Research impact statement

fastDiversity has been used in conservation genomics analyses of threatened Australian flora. @mcmaster2025homoranthus used the package to compare diversity and population structure across multiple species of Homoranthus (Myrtaceae). @mcmaster2025zieria applied the package in a conservation genomic assessment of the critically endangered Zieria obcordata (Rutaceae), identifying extremely low within-population diversity and informing subspecies delimitation. The package has additionally been used in over 30 confidential threatened-species assessments conducted for the New South Wales Government Saving our Species program, demonstrating applicability in operational conservation management contexts.

# Acknowledgements

Development of fastDiversity benefited from discussions with collaborators in conservation genomics and population genetics at the University of Sydney, the Botanic Gardens of Sydney, and Queensland Herbarium.

# AI usage disclosure

Claude Sonnet 4.6 (Anthropic) and ChatGPT (OpenAI GPT-5.5) were used to assist with editing and restructuring documentation and manuscript text based on author-provided drafts and outlines. All scientific content, software design decisions, implementations, analyses, and interpretations were developed, verified, and approved by the author.

# References
