library(testthat)

test_that("individual_F returns a named vector", {
  data(example_data)
  result <- individual_F(example_gt, example_meta$population,
                         maf = 0.05, max_missingness = 0.3)
  expect_true(is.vector(result))
  expect_true(!is.null(names(result)))
})

test_that("individual_F and individual_Ho are negatively correlated", {
  data(example_data)
  ho <- as.numeric(individual_Ho(example_gt, example_meta$population,
                                 maf = 0.05, max_missingness = 0.3))
  f  <- as.numeric(individual_F(example_gt, example_meta$population,
                                maf = 0.05, max_missingness = 0.3))
  expect_true(cor(ho, f, use = "complete.obs") < 0)
})