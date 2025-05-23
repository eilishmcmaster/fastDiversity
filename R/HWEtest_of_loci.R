
HWEtest_of_loci <- function(gt){
  get_allele_freq_matrix <- function(geno_vector) { # make input matrix for HWxtest
    # Fast genotype counting with fixed order: 0=AA, 1=AB, 2=BB
    counts <- tabulate(factor(geno_vector, levels = c(0, 1, 2)), nbins = 3)
    # Assign counts directly: counts = c(AA, AB, BB)
    # Fill the symmetric matrix
    matrix(
      c(counts[1], counts[2],
        NA, counts[3]),
      nrow = 2, dimnames = list(c("A", "B"), c("A", "B"))
    )
  }
  
  # Apply hwx.test per locus
  hwx_results <- lapply(1:ncol(gt), function(index) {
    geno_counts <- get_allele_freq_matrix(gt[[index]])
    result <- tryCatch({
      HWxtest::hwx.test(geno_counts)
    }, error = function(e) {
      warning(paste("Error at locus", index, ":", e$message))
      return(NULL)
    })
    return(result)
  })
  
  names(hwx_results) <- names(gt)
  
  # # get matrix of all p values
  # pvalues_matrix <- do.call(rbind, lapply(hwx_results, function(res) {
  #   if (!is.null(res) && "Pvalues" %in% names(res)) {
  #     res$Pvalues
  #   } else {
  #     rep(NA, length(hwx_results[[which(!sapply(hwx_results, is.null))[1]]]$Pvalues))  # match length
  #   }
  # })) %>% as.data.table()
  
  # Engels 2009, get all LLR P values
  LLR_p <- lapply(hwx_results, function(x) x[["Pvalues"]][["LLR"]]) %>% unlist
  #  U-score test for homozygote or heterozygote excess (Rousset and Raymond 1995), which can be thought as a “one-sided” procedure for narrowing the alternative hypotheses -- very similar to F, where <0 is heterozygote excess, >0 is homozygote excess
  # ExactoHW reports either P(u>observed) or P(u<observed) depending on whether the observed U-score is positive or negative
  U_score <- lapply(hwx_results, function(x) x[["observed"]][["U"]]) %>% unlist
  
  ## Fisher global P value of LLR 
  LLR_X <- -2 * sum(log(LLR_p))
  LLR_df <- 2 * length(LLR_p)
  fisher_combined_p <- pchisq(LLR_X, LLR_df, lower.tail = F)
  
  # return(c('LLR_fisher_p'= fisher_combined_p, 'mean_U_score' = mean(U_score), "loci_HWX_results"=list(hwx_results)))
  return(c('LLR_fisher_p'= fisher_combined_p, 'mean_U_score' = mean(U_score), 
           'loci_HWX_results' = list(data.table(LLR_p, U_score, 'locus' = names(gt)))))
}
