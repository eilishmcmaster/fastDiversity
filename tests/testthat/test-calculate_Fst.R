test_that("calculate_Fst matches hierfstat::fs.dosage population-specific Fst", {
  skip_if_not_installed("hierfstat")
  data(example_data)
  gt <- as.matrix(example_gt)  # includes missing data

  ours <- calculate_Fst(gt, example_meta$population)
  theirs <- hierfstat::fs.dosage(gt, example_meta$population)$Fs["Fst", names(ours)]
  expect_equal(ours, theirs, tolerance = 1e-10)

  keep <- example_meta$population == 1
  ours <- calculate_Fst(gt[keep, ], example_meta$site[keep])
  theirs <- hierfstat::fs.dosage(gt[keep, ], example_meta$site[keep])$Fs["Fst", names(ours)]
  expect_equal(ours, theirs, tolerance = 1e-10)
})

test_that("calculate_Fst is near zero for samples from one population and high for fixed differences", {
  set.seed(1)
  p <- runif(500, 0.1, 0.9)
  one_pop <- t(replicate(60, rbinom(500, 2, p)))
  fst <- calculate_Fst(one_pop, rep(c("A", "B", "C"), each = 20))
  expect_true(all(abs(fst) < 0.02))

  fixed <- rbind(matrix(0, 10, 50), matrix(2, 10, 50))
  expect_equal(unname(calculate_Fst(fixed, rep(c("A", "B"), each = 10))), c(1, 1))
})

test_that("calculate_Fst requires at least two sites", {
  expect_error(calculate_Fst(toy_gt, rep("s", 4)), "at least two sites")
})
