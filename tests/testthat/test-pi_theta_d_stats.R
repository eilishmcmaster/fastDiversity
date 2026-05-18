test_that("pi_theta_d_stats returns a data frame", {
  data(example_data)
  result <- pi_theta_d_stats(
    example_gt,
    genetic_group_variable = example_meta$population,
    site_variable          = example_meta$population,
    minimum_n              = 5,
    minimum_loci           = 5,
    max_missingness        = 0.1,
    maf                    = 0,
    allele_count_min       = 0
  )
  expect_true(is.data.frame(result))
})

test_that("pi_theta_d_stats has expected columns", {
  data(example_data)
  result <- pi_theta_d_stats(
    example_gt,
    genetic_group_variable = example_meta$population,
    site_variable          = example_meta$population,
    minimum_n              = 5,
    minimum_loci           = 5,
    max_missingness        = 0.1,
    maf                    = 0,
    allele_count_min       = 0
  )
  expect_true(all(c("group", "site", "pi", "S", "theta", "D", "loci", "n") %in% colnames(result)))
})

test_that("pi_theta_d_stats numeric columns are numeric", {
  data(example_data)
  result <- pi_theta_d_stats(
    example_gt,
    genetic_group_variable = example_meta$population,
    site_variable          = example_meta$population,
    minimum_n              = 5,
    minimum_loci           = 5,
    max_missingness        = 0.1,
    maf                    = 0,
    allele_count_min       = 0
  )
  expect_true(all(!is.na(as.numeric(result$pi))))
  expect_true(all(!is.na(as.numeric(result$theta))))
  expect_true(all(as.numeric(result$pi) >= 0))
  expect_true(all(as.numeric(result$theta) >= 0))
  expect_true(all(as.numeric(result$S) >= 0))
})

test_that("pi_theta_d_stats returns one row per group when site equals group", {
  data(example_data)
  result <- pi_theta_d_stats(
    example_gt,
    genetic_group_variable = example_meta$population,
    site_variable          = example_meta$population,
    minimum_n              = 5,
    minimum_loci           = 5,
    max_missingness        = 0.1,
    maf                    = 0,
    allele_count_min       = 0
  )
  expect_equal(nrow(result), length(unique(example_meta$population)))
})

test_that("pi_theta_d_stats n argument subsamples correctly", {
  data(example_data)
  result <- pi_theta_d_stats(
    example_gt,
    genetic_group_variable = example_meta$population,
    site_variable          = example_meta$population,
    n                      = 10,
    minimum_n              = 5,
    minimum_loci           = 5,
    max_missingness        = 0.1,
    maf                    = 0,
    allele_count_min       = 0
  )
  expect_true(all(as.numeric(result$n) <= 10))
})

test_that("pi_theta_d_stats returns NULL when minimum_n is too high", {
  data(example_data)
  expect_warning(
    result <- pi_theta_d_stats(
      example_gt,
      genetic_group_variable = example_meta$population,
      site_variable          = example_meta$population,
      minimum_n              = 9999,
      minimum_loci           = 5,
      max_missingness        = 0.1,
      maf                    = 0,
      allele_count_min       = 0
    ),
    "No data passed filters"
  )
  expect_null(result)
})

test_that("pi_theta_d_stats S increases with sample size", {
  data(example_data)
  result_small <- pi_theta_d_stats(
    example_gt,
    genetic_group_variable = example_meta$population,
    site_variable          = example_meta$population,
    n                      = 5,
    minimum_n              = 5,
    minimum_loci           = 5,
    max_missingness        = 0.1,
    maf                    = 0,
    allele_count_min       = 0
  )
  result_large <- pi_theta_d_stats(
    example_gt,
    genetic_group_variable = example_meta$population,
    site_variable          = example_meta$population,
    n                      = NULL,
    minimum_n              = 5,
    minimum_loci           = 5,
    max_missingness        = 0.1,
    maf                    = 0,
    allele_count_min       = 0
  )
  # S should be equal or greater with more individuals on average
  expect_true(mean(as.numeric(result_large$S)) >= mean(as.numeric(result_small$S)))
})