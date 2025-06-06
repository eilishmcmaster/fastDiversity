% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/individual_Ho.R
\name{individual_Ho}
\alias{individual_Ho}
\title{Calculate individual heterozygosity}
\usage{
individual_Ho(
  gt,
  genetic_group_variable = NULL,
  minimum_loci = 50,
  maf = 0.05,
  max_missingness = 0.3
)
}
\arguments{
\item{gt}{Genotype matrix where individuals are rows and loci are columns, coded as 
0=aa, 1=aA, 2=AA.}

\item{genetic_group_variable}{(optional) Vector of genetic group IDs for each individual (same 
length as number of rows in gt). Can represent species, 
subspecies, or other biologically relevant structure groups. If not
included, the whole gt will be treated as one group.}

\item{minimum_loci}{(optional) Minimum number of loci required to proceed with analysis for each 
genetic group (genetic groups with fewer loci will be skipped).}

\item{maf}{(optional) Minimum minor allele frequency filter.}

\item{max_missingness}{(optional) Maximum proportion of missing data allowed per locus (loci with 
higher missingness will be filtered out).}
}
\value{
A named vector containing F for each individual in gt.
}
\description{
Calculate individual heterozygosity
}
\examples{
data(example_data)
# get inbreeding coefficients for individuals by genetic group
Hos <- individual_Ho(example_gt, example_meta$population, maf=0.05)
}
\references{
Frankham, R., Ballou, J. D., & Briscoe, D. A. (2010). Introduction to Conservation Genetics (2nd ed.). Cambridge University Press. https://doi.org/10.1017/CBO9780511809002
}
