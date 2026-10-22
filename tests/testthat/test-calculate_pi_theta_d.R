test_that("calculate_pi_theta_d matches hand-calculated pi and Watterson's theta", {
  res <- calculate_pi_theta_d(toy_gt)

  # pi per locus = 2pq * n/(n-1), with n = number of sampled alleles (8, 6, 8)
  expected_pi <- (2 * 1/2 * 1/2 * 8/7) + (2 * 1/6 * 5/6 * 6/5) + (2 * 7/8 * 1/8 * 8/7)
  # theta = sum over segregating sites of 1 / a_n, a_n = sum_{i=1}^{n-1} 1/i
  a_n <- function(n) sum(1 / seq_len(n - 1))
  expected_theta <- 1 / a_n(8) + 1 / a_n(6) + 1 / a_n(8)

  expect_equal(res$raw_pi, expected_pi)
  expect_equal(res$watterson_theta, expected_theta)
  expect_equal(res$num_sites, 3)
})

test_that("calculate_pi_theta_d ignores monomorphic loci", {
  with_mono <- cbind(toy_gt, M1 = c(0, 0, 0, 0), M2 = c(2, 2, NA, 2))
  expect_equal(calculate_pi_theta_d(with_mono), calculate_pi_theta_d(toy_gt))
})
