library(testthat)

test_that("make_allele_list returns a named list", {
  data(example_data)
  result <- make_allele_list(example_gt, example_meta$population)
  expect_true(is.list(result))
  expect_true(!is.null(names(result)))
})

test_that("make_allele_list has one entry per group", {
  data(example_data)
  result <- make_allele_list(example_gt, example_meta$population)
  expect_equal(length(result), length(unique(example_meta$population)))
})