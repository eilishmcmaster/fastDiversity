#' Create Allele List by Group
#' 
#' This function creates a list of alleles for each group specified by `site_variable` 
#' based on genotype data (`gt`). It iterates over unique groups, calculates alleles 
#' for each group, and returns a list where each element corresponds to a group.
#' 
#' @param gt Genotype matrix where rows represent individuals and columns represent loci (0, 1, 2 format).
#' @param site_variable Vector specifying group IDs (sites) for each individual.
#' @return List where each element contains a vector of alleles for a group.
#' @export
make_allele_list <- function(gt, site_variable) {
  # Get unique groups
  groups <- unique(site_variable)
  
  # Initialize output list
  out <- vector("list", length(groups))
  names(out) <- groups
  
  # Prepare loci data frame outside the loop
  loci <- data.frame("loci" = colnames(gt),
                     "allele1" = paste0(1:ncol(gt), "_a"),
                     "allele2" = paste0(1:ncol(gt), "_b"))
  
  # Iterate over each group and calculate alleles
  for (i in seq_along(groups)) {
    group_mask <- site_variable == groups[i]
    df <- gt[group_mask, , drop = FALSE] # Subset without dropping dimensions
    cat("There are ", nrow(df), " samples in ", groups[i], "\n")
    df2 <- rbind("names" = colnames(df), df)
    alleles <- apply(df2, 2, matcher, loci)
    out[[i]] <- unlist(alleles, use.names = FALSE)
  }
  
  return(out)
}
