#' Total Allele Counts
#' 
#' This function calculates the number of total and private alleles unique to each group based on 
#' biallelic gentoype data. It computes the count of private alleles and total alleles 
#' for each group and aggregates results into a data table.
#' 
#' @param gt Genotype matrix where individuals are rows and loci are columns, coded as 
#'           0=aa, 1=aA, 2=AA.
#' @param group_variable Vector of genetic group IDs for each individual (same 
#'                               length as number of rows in gt). Can represent sites, populations etc.
#' @param maf Minimum minor allele frequency filter.
#' @return Data table summarizing private alleles and total alleles for each group and a total count.
#' @export
calculate_total_alleles <- function(gt, group_variable, maf) {
  groups <- unique({{group_variable}})
  out <- vector("list", length(groups))  # Initialize output as a list with the same length as groups
  
  ds <- gt
  keepers <- get_minor_allele_frequencies(ds)
  ds <- ds[, which(keepers >= maf)]
  cat("Found ", ncol(ds), " poly sites\n")     
  
  ds_inverse <- 2 - ds
  colnames(ds) <- paste(colnames(ds), 'A', sep="_")
  colnames(ds_inverse) <- paste(colnames(ds_inverse), 'B', sep="_")
  x <- cbind(ds, ds_inverse)
  
  x2 <- x %>% as.matrix()
  x2[x2==2] <- 1
  n <- c()
  out <- list()
  
  for (i in seq_along(groups)) {
    # ds[{{group_variable}} %in% groups[i], , drop = FALSE]
    x3 <- x2[{{group_variable}} %in% groups[i], , drop = FALSE]
    cat("There are ", nrow(x3), " samples in ", groups[i], "\n")
    # x3 <- rbind("names" = colnames(x3), as.data.frame(x3))
    
    n <- c(n, nrow(x3))
    allele_max <- colSums(x3, na.rm=TRUE) %>% as.vector()
    allele_max[allele_max>=1] <- 1 # presence absence
    out[[i]] <- as.vector(allele_max)
    names(out)[i] <- groups[i]
  }
  
  z <- do.call(rbind,out)
  total_allele_count <- rowSums(z, na.rm=TRUE)
  pa_loci <- which(colSums(z, na.rm=TRUE)==1)
  
  if(length(pa_loci)>1){
    private_allele_count <- rowSums(z[,pa_loci], na.rm=TRUE)
  }
  if(length(pa_loci)==1){
    private_allele_count <- rep(0,nrow(z))
    private_allele_count[which(z[,pa_loci]==1)] <- 1
  }
  if(length(pa_loci)==0){
    private_allele_count <- rep(0,nrow(z))
  }
  out_df <- cbind(private_allele_count, total_allele_count, n) %>% as.data.frame()
  out_df$population <- rownames(out_df)
  rownames(out_df) <- NULL
  out_df <- rbind(out_df, c(ncol(x2), ncol(x2),sum(n), 'Total'))
  out_df <- out_df[,c(4,1,2,3)]
  return(out_df)
}
