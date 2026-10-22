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

test_that("faststats returns numeric statistic columns", {
  data(example_data)
  invisible(capture.output(
    result <- faststats(example_gt, example_meta$population, example_meta$site)
  ))
  stat_cols <- setdiff(names(result), c("group", "site"))
  expect_true(all(vapply(result[stat_cols], is.numeric, logical(1))))
})

test_that("faststats matches hand-calculated values for a single site", {
  invisible(capture.output(
    result <- faststats(toy_gt, rep("g", 4), rep("s", 4),
                        minimum_n = 3, minimum_loci = 1, maf = 0, max_missingness = 1)
  ))
  He <- mean(toy_Hes)
  expect_equal(result$Ho, round(4/11, 4))
  expect_equal(result$He, round(He, 4))
  expect_equal(result$uHe, round(mean(2 * toy_ns / (2 * toy_ns - 1) * toy_Hes), 4))
  expect_equal(result$Fis, round(1 - round(4/11, 4) / round(He, 4), 4))
  expect_equal(result$Fst, 0)  # one site, so Hs equals Ht
  expect_equal(result$loci, 3)
  expect_equal(result$n, 4)
})

test_that("faststats return_locus_stats works without run_HWE_test", {
  data(example_data)
  keep <- example_meta$population == 2
  expect_message(
    invisible(capture.output(
      result <- faststats(example_gt[keep, ], example_meta$population[keep],
                          example_meta$site[keep], return_locus_stats = TRUE)
    )),
    "setting run_HWE_test = TRUE"
  )
  expect_named(result, c("summary", "locus_HWE"))
  expect_true("HWE_p" %in% names(result$summary))
})

test_that("faststats adds bootstrap CI columns when requested", {
  data(example_data)
  keep <- example_meta$population == 2
  set.seed(1)
  invisible(capture.output(
    result <- faststats(example_gt[keep, ], example_meta$population[keep],
                        example_meta$site[keep], boots = 10,
                        get_CI_loci_resampling = TRUE,
                        get_CI_individual_resampling = TRUE)
  ))
  expect_true(all(c("loci_ho.2.5%", "loci_ho.97.5%",
                    "individuals_ho.2.5%", "individuals_ho.97.5%") %in% names(result)))
  expect_true(all(result$`loci_ho.2.5%` <= result$`loci_ho.97.5%`))
})
