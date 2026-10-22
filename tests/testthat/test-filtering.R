# Hand-worked example combining all faststats() filters.
#
# Group g1 has 8 individuals at three sites; group g2 has 2 individuals.
#
#          L1  L2  L3  L4  L5
#   s1a     0   0   0   0   0
#   s1b     1   0  NA   0   1
#   s1c     2   0   1   0   1
#   s2a     0   0  NA   0   2
#   s2b     0   0  NA   1   2
#   s2c     1   0   0   0  NA
#   s3a     2   0   1   0   0
#   s3b     1   0  NA   0   1
#
# Filters with minimum_n = 3, maf = 0.1, max_missingness = 0.3, applied within g1:
#   g2 is dropped (2 individuals < minimum_n)
#   s3 is dropped (2 individuals < minimum_n), but its individuals still count
#     towards locus filtering and Ht for the group
#   L2 is dropped (monomorphic, MAF = 0)
#   L3 is dropped (4/8 = 50% missing)
#   L4 is dropped (MAF = 1/16 = 0.0625)
#   L1 (alt freq 7/16) and L5 (alt freq 7/14, 1/8 missing) are kept
filter_gt <- matrix(
  c(0, 1, 2, 0, 0, 1, 2, 1,   0, 1,
    0, 0, 0, 0, 0, 0, 0, 0,   1, 1,
    0, NA, 1, NA, NA, 0, 1, NA, 0, 1,
    0, 0, 0, 0, 1, 0, 0, 0,   0, 1,
    0, 1, 1, 2, 2, NA, 0, 1,  2, 2),
  nrow = 10,
  dimnames = list(c("s1a", "s1b", "s1c", "s2a", "s2b", "s2c", "s3a", "s3b", "g2a", "g2b"),
                  paste0("L", 1:5))
)
filter_groups <- c(rep("g1", 8), "g2", "g2")
filter_sites  <- c(rep("s1", 3), rep("s2", 3), rep("s3", 2), "s4", "s4")

run_filter_example <- function(minimum_loci = 2) {
  suppressMessages(faststats(filter_gt, filter_groups, filter_sites,
                             minimum_n = 3, minimum_loci = minimum_loci,
                             maf = 0.1, max_missingness = 0.3))
}

test_that("faststats drops small groups and sites and keeps only loci passing filters", {
  result <- run_filter_example()
  expect_equal(result$group, c("g1", "g1"))
  expect_equal(result$site, c("s1", "s2"))
  expect_equal(result$loci, c(2, 2))  # L1 and L5
  expect_equal(result$n, c(3, 3))
})

test_that("faststats statistics after filtering match hand calculations", {
  result <- run_filter_example()

  # s1 on L1 = (0, 1, 2) and L5 = (0, 1, 1)
  #   Ho = 3 hets / 6 genotypes; He = mean(1/2, 4/9); uHe = mean(6/5 * c(1/2, 4/9))
  # s2 on L1 = (0, 0, 1) and L5 = (2, 2, NA)
  #   Ho = 1 het / 5 genotypes; He = mean(10/36, 0); uHe = mean(6/5 * 10/36, 4/3 * 0)
  expect_equal(result$Ho,  round(c(3/6, 1/5), 4))
  expect_equal(result$He,  round(c(mean(c(1/2, 4/9)), mean(c(10/36, 0))), 4))
  expect_equal(result$uHe, round(c(mean(6/5 * c(1/2, 4/9)), mean(c(6/5 * 10/36, 0))), 4))

  # Ht uses all 8 g1 individuals, including s3: L1 p = 7/16, L5 p = 1/2
  Ht <- mean(c(2 * 7/16 * 9/16, 2 * 1/2 * 1/2))
  He <- c(mean(c(1/2, 4/9)), mean(c(10/36, 0)))
  expect_equal(result$Gst, round((Ht - He) / Ht, 4))

  # Rarefied Ar across s1 and s2: L1 rarefied to 6 alleles, L5 to 4 alleles
  #   s1: L1 = 2, L5 = 1 + (1 - choose(4, 4) / choose(6, 4));  s2: L1 = 2, L5 = 1
  expect_equal(result$Ar, round(c(mean(c(2, 1 + (1 - 1 / choose(6, 4)))), mean(c(2, 1))), 4))
})

test_that("faststats skips a group when too few loci pass the filters", {
  expect_message(
    result <- faststats(filter_gt, filter_groups, filter_sites, minimum_n = 3,
                        minimum_loci = 3, maf = 0.1, max_missingness = 0.3),
    "does not have enough loci"
  )
  expect_equal(nrow(result), 0)
})

test_that("relaxing each filter retains the expected extra locus", {
  base <- run_filter_example()
  no_maf <- suppressMessages(faststats(filter_gt, filter_groups, filter_sites, minimum_n = 3,
                                       minimum_loci = 2, maf = 0.05, max_missingness = 0.3))
  no_missing <- suppressMessages(faststats(filter_gt, filter_groups, filter_sites, minimum_n = 3,
                                           minimum_loci = 2, maf = 0.1, max_missingness = 0.5))
  expect_equal(no_maf$loci, base$loci + 1)      # L4 (MAF 0.0625) now passes
  expect_equal(no_missing$loci, base$loci + 1)  # L3 (50% missing) now passes
})
