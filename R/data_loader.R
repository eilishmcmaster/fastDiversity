#' @title Load Sample Dataset 1
#' @description This function loads the first sample dataset stored in CSV format.
#' @name example_gt
#' @docType data
#' @usage data(example_gt)
#' @examples
#' data(example_gt)
#' head(example_gt)
"example_gt" <- read.csv(system.file("data", "example_gt.csv", package = "fastDiversity"), row.names = 1)

#' @title Load Sample Dataset 2
#' @description This function loads the second sample dataset stored in CSV format.
#' @name example_meta
#' @docType data
#' @usage data(example_meta)
#' @examples
#' data(example_meta)
#' head(example_meta)
"example_meta" <- read.csv(system.file("data", "example_meta.csv", package = "fastDiversity"))
