library(testthat)

test_that("individual_Ho returns a named vector", {
  data(example_data)
  result <- individual_Ho(example_gt, example_meta$population,
                          maf = 0.05, max_missingness = 0.3)
  expect_true(is.vector(result))
  expect_true(!is.null(names(result)))
})

test_that("individual_Ho values are between 0 and 1", {
  data(example_data)
  result <- individual_Ho(example_gt, example_meta$population,
                          maf = 0.05, max_missingness = 0.3)
  expect_true(all(as.numeric(result) >= 0 & as.numeric(result) <= 1, na.rm = TRUE))
})

test_that("individual_Ho returns one value per sample", {
  data(example_data)
  result <- individual_Ho(example_gt, example_meta$population,
                          maf = 0.05, max_missingness = 0.3)
  expect_equal(length(result), nrow(example_gt))
})