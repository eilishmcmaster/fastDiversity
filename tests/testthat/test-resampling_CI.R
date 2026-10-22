test_that("resampling_loci_CI returns ordered intervals for F, He and Ho", {
  data(example_data)
  set.seed(1)
  ci <- resampling_loci_CI(example_gt[1:30, ], boots = 50)

  expect_named(ci, c(paste0("loci_f.", c("2.5%", "97.5%", "mean")),
                     paste0("loci_he.", c("2.5%", "97.5%", "mean")),
                     paste0("loci_ho.", c("2.5%", "97.5%", "mean"))))
  for (stat in c("loci_f", "loci_he", "loci_ho")) {
    lower <- ci[[paste0(stat, ".2.5%")]]
    upper <- ci[[paste0(stat, ".97.5%")]]
    mid   <- ci[[paste0(stat, ".mean")]]
    expect_lte(lower, mid)
    expect_lte(mid, upper)
  }
})

test_that("resampling_individuals_CI returns ordered intervals for F, He and Ho", {
  data(example_data)
  set.seed(1)
  ci <- resampling_individuals_CI(example_gt[1:30, ], boots = 50)

  for (stat in c("individuals_f", "individuals_he", "individuals_ho")) {
    lower <- ci[[paste0(stat, ".2.5%")]]
    upper <- ci[[paste0(stat, ".97.5%")]]
    mid   <- ci[[paste0(stat, ".mean")]]
    expect_lte(lower, mid)
    expect_lte(mid, upper)
  }
})

test_that("bootstrap He and Ho intervals contain the point estimates", {
  data(example_data)
  gt <- example_gt[1:30, ]
  set.seed(1)
  ci <- resampling_loci_CI(gt, boots = 200)
  # mean of per-locus values, as calculated in the bootstrap
  ho_hat <- mean(colMeans(gt == 1, na.rm = TRUE), na.rm = TRUE)
  he_hat <- mean(calculate_Hes(gt), na.rm = TRUE)

  expect_gte(he_hat, ci[["loci_he.2.5%"]])
  expect_lte(he_hat, ci[["loci_he.97.5%"]])
  expect_gte(ho_hat, ci[["loci_ho.2.5%"]])
  expect_lte(ho_hat, ci[["loci_ho.97.5%"]])
})

test_that("resampling is reproducible with set.seed and does not reset the user's seed", {
  data(example_data)
  gt <- example_gt[1:30, 1:200]

  set.seed(42); a <- resampling_loci_CI(gt, boots = 20)
  set.seed(42); b <- resampling_loci_CI(gt, boots = 20)
  expect_identical(a, b)

  # consecutive calls should draw different resamples (a fixed internal seed made them identical)
  set.seed(42)
  first  <- resampling_loci_CI(gt, boots = 20)
  second <- resampling_loci_CI(gt, boots = 20)
  expect_false(identical(first, second))

  first  <- resampling_individuals_CI(gt, boots = 20)
  second <- resampling_individuals_CI(gt, boots = 20)
  expect_false(identical(first, second))
})

test_that("resampling a matrix with no variation gives zero-width intervals", {
  gt <- matrix(rep(c(0, 1, 2, 1), 10), nrow = 4)  # identical loci
  set.seed(1)
  ci <- resampling_loci_CI(gt, boots = 20)
  expect_equal(ci[["loci_ho.2.5%"]], ci[["loci_ho.97.5%"]])
  expect_equal(ci[["loci_he.2.5%"]], ci[["loci_he.97.5%"]])
})
