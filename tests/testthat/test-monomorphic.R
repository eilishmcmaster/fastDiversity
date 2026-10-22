# How monomorphic (invariant) sites affect each statistic.

test_that("pi, theta and Tajima's D do not change when monomorphic sites are added", {
  data(example_data)
  gt <- as.matrix(example_gt[example_meta$population == 1, ])
  invariant <- matrix(rep(c(0, 2), each = nrow(gt), times = 250), nrow(gt))  # all-0 and all-2 columns
  colnames(invariant) <- paste0("m", 1:500)

  snps_only <- calculate_pi_theta_d(gt)
  with_invariant <- calculate_pi_theta_d(cbind(gt, invariant))

  expect_equal(with_invariant, snps_only)
})

test_that("faststats ratio statistics do not depend on monomorphic loci, but Ho and He do", {
  data(example_data)
  keep <- example_meta$population == 1
  gt <- as.matrix(example_gt[keep, ])
  gt <- gt[, colSums(is.na(gt)) == 0]  # complete data
  sites <- example_meta$site[keep]
  invariant <- matrix(0, nrow(gt), 1000, dimnames = list(NULL, paste0("m", 1:1000)))

  run <- function(g) suppressMessages(faststats(g, rep("1", nrow(g)), sites, maf = 0,
                                                max_missingness = 1, minimum_loci = 1))
  snps_only <- run(gt)
  with_invariant <- run(cbind(gt, invariant))

  for (stat in c("Fis", "uFis", "Gst", "Fst")) {
    # equal up to rounding of the output to 4 decimal places
    expect_true(all(abs(with_invariant[[stat]] - snps_only[[stat]]) <= 1e-4 + 1e-9, na.rm = TRUE))
  }
  # Ho and He are averages over loci, so they shrink in proportion to the share of SNPs
  # (only sites with polymorphic loci are compared, since the ratio is undefined otherwise)
  shrink <- ncol(gt) / (ncol(gt) + 1000)
  expect_equal(with_invariant$Ho, snps_only$Ho * shrink, tolerance = 1e-3)
  expect_equal(with_invariant$He, snps_only$He * shrink, tolerance = 1e-3)
})
