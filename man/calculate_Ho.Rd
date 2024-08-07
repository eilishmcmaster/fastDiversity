% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_Ho.R
\name{calculate_Ho}
\alias{calculate_Ho}
\title{Calculate Observed Heterozygosity (Ho)}
\usage{
calculate_Ho(gt)
}
\arguments{
\item{gt}{Genotype matrix where individuals are rows and loci are columns, coded as 
0=aa, 1=aA, 2=AA.}
}
\value{
Observed heterozygosity (Ho) across all loci.
}
\description{
This function calculates observed heterozygosity (Ho) based on genotype data.
Heterozygosity (Ho) is defined as the proportion of heterozygous individuals 
among all individuals with genotype data.
}
