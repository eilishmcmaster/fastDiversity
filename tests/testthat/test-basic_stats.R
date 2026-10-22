test_that("calculate_Ho pools heterozygotes over all non-missing genotypes", {
  # 4 heterozygotes among 11 non-missing genotypes
  expect_equal(calculate_Ho(toy_gt), 4/11)
})

test_that("calculate_Hes returns 2pq per locus", {
  expect_equal(calculate_Hes(toy_gt), toy_Hes)
})

test_that("calculate_uHe applies the 2n/(2n-1) correction per locus", {
  expected <- mean(2 * toy_ns / (2 * toy_ns - 1) * toy_Hes)
  expect_equal(calculate_uHe(toy_ns, toy_Hes), expected)
})

test_that("minor allele frequencies and counts are correct", {
  expect_equal(unname(get_minor_allele_frequencies(toy_gt)), c(1/2, 1/6, 1/8))
  expect_equal(unname(get_minor_allele_counts(toy_gt)), c(4, 1, 1))
})

test_that("individual_Ho returns per-individual heterozygosity", {
  ho <- individual_Ho(toy_gt, minimum_loci = 1, maf = 0, max_missingness = 1)
  expect_equal(ho, c(i1 = 0, i2 = 1/3, i3 = 1/3, i4 = 1))
})

test_that("individual_F is 1 - individual Ho / mean He", {
  f <- individual_F(toy_gt, minimum_loci = 1, maf = 0, max_missingness = 1)
  expected <- 1 - c(0, 1/3, 1/3, 1) / mean(toy_Hes)
  expect_equal(unname(f), expected)
  expect_equal(names(f), rownames(toy_gt))
})

test_that("individual_Ho returns NA instead of erroring when too few loci remain", {
  expect_message(
    ho <- individual_Ho(toy_gt, minimum_loci = 10, maf = 0, max_missingness = 1),
    "not have enough loci"
  )
  expect_true(all(is.na(ho)))
  expect_length(ho, nrow(toy_gt))
})
