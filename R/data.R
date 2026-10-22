#' Example SNP genotype data and sample metadata
#'
#' DArTseq SNP genotypes for 257 individuals of an unpublished plant species
#' collected in eastern Australia, with matching sample metadata. The species
#' identity is withheld pending publication. \code{data(example_data)} loads
#' two objects, \code{example_gt} and \code{example_meta}.
#'
#' @format
#' \code{example_gt}: a data frame with 257 rows (individuals) and 1906 columns
#' (SNP loci), with genotypes coded as 0 (homozygous reference, aa),
#' 1 (heterozygous, aA), 2 (homozygous alternate, AA) and \code{NA} (missing;
#' about 13\% of genotypes). Row names are sample IDs and column names are
#' DArTseq locus IDs. All loci are polymorphic across the full dataset;
#' monomorphic sites are not included.
#'
#' \code{example_meta}: a data frame with 257 rows, in the same order as the
#' rows of \code{example_gt}, and 3 columns:
#' \describe{
#'   \item{sample}{Sample ID, matching \code{rownames(example_gt)}.}
#'   \item{site}{Collection site (55 sites).}
#'   \item{population}{Genetic group assignment (integer, 5 groups).}
#' }
#' @source Unpublished DArTseq data provided by the package author.
#' @docType data
#' @name example_data
#' @aliases example_gt example_meta
#' @usage data(example_data)
#' @keywords datasets
#' @examples
#' data(example_data)
#' dim(example_gt)
#' example_gt[1:5, 1:5]
#' head(example_meta)
#' table(example_meta$population)
NULL
