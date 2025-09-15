#' This package provides functions for basic population genetics analysis using SNP data
#' in biallelic format (0, 1, 2). It includes methods for calculating diversity statistics,
#' allelic richness, heterozygosity, and more.
#'
#' @title fastDiversity: Population Genetics Analysis
#' @description
#' This package was developed to simplify population genetic analysis tasks by offering
#' efficient functions compatible with SNP data in biallelic format (0, 1, 2), especially
#' when there are multiple species or subpopulations in a single dataset.
#'
#' @details
#' The package includes functions for:
#' - Calculating basic diversity statistics.
#' - Estimating allelic richness using rarefaction.
#' - Identifying private alleles across populations.
#'
#' For more details on each function, use `?function_name`.
#' 
#' @docType _PACKAGE
#' @name fastDiversity
#' @seealso \code{\link{faststats}}, \code{\link{calculate_private_alleles}}
#' @references
#' Frankham, R., Ballou, J. D., & Briscoe, D. A. (2010). Introduction to Conservation Genetics (2nd ed.). Cambridge University Press. https://doi.org/10.1017/CBO9780511809002
#' 
#' Kalinowski ST (2005). "HP-Rare." Molecular Ecology Notes. URL: https://www.montana.edu/kalinowski/software/documents/2005_HP-Rare_MolecularEcologyNotes.pdf
#' 
#' Keenan K, McGinnity P, Cross TF, et al. (2013). "diveRsity: An R package for the estimation of population genetics parameters and their associated errors." Methods in Ecology and Evolution. URL: https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210x.12067
#'
#' @author
#' Eilish S. McMaster <emcm4052@uni.sydney.edu.au> <eilish.mcmaster@botanicgardens.nsw.gov.au>
#' @keywords population-genetics package
NULL