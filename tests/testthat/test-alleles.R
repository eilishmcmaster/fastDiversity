test_that("make_allele_list assigns alleles using 0=aa, 1=aA, 2=AA coding", {
  gt <- matrix(c(0, 0,    # L1: aa only      -> allele a
                 1, 1,    # L2: aA only      -> alleles a and b
                 2, 2),   # L3: AA only      -> allele b
               nrow = 2, dimnames = list(c("i1", "i2"), c("L1", "L2", "L3")))
  invisible(capture.output(res <- make_allele_list(gt, c("pop", "pop"))))
  expect_setequal(res$pop, c("1_a", "2_a", "2_b", "3_b"))
})

test_that("calculate_private_alleles counts private and total alleles per group", {
  allele_list <- list(A = c("1_a", "1_b", "2_a"),
                      B = c("1_a", "2_a", "2_b", "3_a"))
  res <- calculate_private_alleles(allele_list)

  expect_equal(res$population, c("A", "B", "Total"))
  expect_equal(as.numeric(res$private_allele_count), c(1, 2, 5))  # A: 1_b; B: 2_b, 3_a
  expect_equal(as.numeric(res$total_allele_count),   c(3, 4, 5))
})

test_that("make_allele_list and calculate_total_alleles agree on private and total alleles", {
  data(example_data)
  invisible(capture.output({
    a <- calculate_private_alleles(make_allele_list(example_gt, example_meta$population))
    b <- calculate_total_alleles(example_gt, example_meta$population, maf = 0)
  }))
  expect_equal(as.numeric(a$private_allele_count), as.numeric(b$private_allele_count))
  expect_equal(as.numeric(a$total_allele_count),   as.numeric(b$total_allele_count))
})

test_that("calculate_rarified_Ar matches hand-calculated rarefaction", {
  # site S1 has 2 individuals, S2 has 3, so rarefy to g = 4 alleles per locus
  gt <- matrix(c(0, 1, 0, 0, 0,     # L1: S1 has 1 alt / 3 ref, S2 has 0 alt / 6 ref
                 2, 2, 1, 1, 2),    # L2: S1 has 4 alt / 0 ref, S2 has 4 alt / 2 ref
               nrow = 5)
  sites <- c("S1", "S1", "S2", "S2", "S2")
  ar <- calculate_rarified_Ar(gt, sites)

  # P(allele seen in 4 draws without replacement) = 1 - choose(N - k, 4) / choose(N, 4)
  # L1, S1: both alleles certain (2); S2: only ref present (1)
  # L2, S1: only alt present (1);   S2: alt certain, ref seen with prob 1 - choose(4,4)/choose(6,4)
  expected <- matrix(c(2, 1,
                       1, 1 + (1 - 1 / choose(6, 4))),
                     nrow = 2)
  expect_equal(ar, expected)
})
