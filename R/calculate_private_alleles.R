#' Calculate Private Alleles
#' 
#' This function calculates the number of private alleles unique to each group based on 
#' allele lists (`allele_list`). It computes the count of private alleles and total alleles 
#' for each group and aggregates results into a data table.
#' 
#' @param allele_list List where each element contains a vector of alleles for a group.
#' @return Data table summarizing private alleles and total alleles for each group and a total count.
#' @export
calculate_private_alleles <- function(allele_list) {
  # Combine all alleles into a single vector
  all_alleles <- unlist(allele_list)
  
  # Calculate total alleles and unique alleles
  total_alleles <- unique(all_alleles)
  result_table <- data.frame(population = character(length(allele_list) + 1), 
                             private_allele_count = numeric(length(allele_list) + 1), 
                             total_allele_count = numeric(length(allele_list) + 1))
  
  for (i in seq_along(allele_list)) {
    current_pop <- unique(allele_list[[i]])
    other_pops <- unlist(allele_list[-i])
    private_alleles <- setdiff(current_pop, other_pops)
    
    result_table[i, ] <- data.frame(population = names(allele_list)[i], 
                                    private_allele_count = length(private_alleles), 
                                    total_allele_count = length(current_pop))
  }
  
  result_table[nrow(result_table), ] <- data.frame(population = "Total", 
                                                   private_allele_count = length(total_alleles), 
                                                   total_allele_count = length(total_alleles))
  
  return(result_table)
}
