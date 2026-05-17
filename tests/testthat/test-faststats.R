library(testthat)

test_that("faststats returns a data frame", {
  data(example_data)
  result <- faststats(example_gt, example_meta$population,
                      example_meta$site, minimum_n = 3,
                      minimum_loci = 50, maf = 0.05, max_missingness = 0.3)
  expect_true(is.data.frame(result))
})

test_that("faststats has expected columns", {
  data(example_data)
  result <- faststats(example_gt, example_meta$population,
                      example_meta$site, minimum_n = 3,
                      minimum_loci = 50, maf = 0.05, max_missingness = 0.3)
  expect_true(all(c("group", "site", "Ho", "He", "Fis", "Fst") %in% colnames(result)))
})

test_that("faststats Ho values are between 0 and 1", {
  data(example_data)
  result <- faststats(example_gt, example_meta$population,
                      example_meta$site, minimum_n = 3,
                      minimum_loci = 50, maf = 0.05, max_missingness = 0.3)
  ho <- as.numeric(result$Ho)
  expect_true(all(ho >= 0 & ho <= 1, na.rm = TRUE))
})

test_that("faststats respects minimum_n filtering", {
  data(example_data)
  # setting minimum_n very high should return fewer or no rows
  result_strict <- faststats(example_gt, example_meta$population,
                             example_meta$site, minimum_n = 999,
                             minimum_loci = 50, maf = 0.05, max_missingness = 0.3)
  expect_true(nrow(result_strict) == 0)
})
