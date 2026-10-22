test_that("MC_HWEtest_of_loci gives a high p-value for a locus in HWE and a low one for heterozygote excess", {
  # 25 AA : 50 AB : 25 BB is exactly at HWE; 0 : 100 : 0 is an extreme heterozygote excess
  hwe_locus    <- c(rep(0, 25), rep(1, 50), rep(2, 25))
  excess_locus <- rep(1, 100)
  res <- MC_HWEtest_of_loci(data.frame(hwe = hwe_locus, excess = excess_locus))

  expect_gt(res$loci_HWX_results[["hwe"]], 0.5)
  expect_lt(res$loci_HWX_results[["excess"]], 1e-6)
})

test_that("MC_HWEtest_of_loci combines locus p-values with Fisher's method", {
  data(example_data)
  gt <- example_gt[1:40, 1:30]
  res <- MC_HWEtest_of_loci(gt)

  p <- res$loci_HWX_results
  expected <- pchisq(-2 * sum(log(p), na.rm = TRUE), 2 * sum(!is.na(p)), lower.tail = FALSE)
  expect_equal(res$HW_fisher_p, round(expected, 6))
})

test_that("MC_HWEtest_of_loci gives the same result for matrix and data frame input", {
  data(example_data)
  gt <- example_gt[1:40, 1:30]
  expect_equal(MC_HWEtest_of_loci(as.matrix(gt)), MC_HWEtest_of_loci(gt))
})

test_that("MC_HWEtest_of_loci does not warn about unlabelled genotypes", {
  expect_no_warning(MC_HWEtest_of_loci(as.data.frame(toy_gt)))
})
