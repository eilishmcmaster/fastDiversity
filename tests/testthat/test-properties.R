# Property-based tests: simulate random dosage matrices, groupings and filter settings,
# and check invariants that should hold for any input.

random_gt <- function(n_ind, n_loci, missing_rate) {
  gt <- vapply(runif(n_loci), function(p) rbinom(n_ind, 2, p), integer(n_ind))
  gt <- matrix(gt, nrow = n_ind)
  gt[runif(length(gt)) < missing_rate] <- NA
  dimnames(gt) <- list(paste0("ind", seq_len(n_ind)), paste0("locus", seq_len(n_loci)))
  gt
}

random_case <- function() {
  n_ind  <- sample(4:80, 1)
  n_loci <- sample(2:400, 1)
  list(
    gt              = random_gt(n_ind, n_loci, missing_rate = runif(1, 0, 0.4)),
    groups          = sample(c("A", "B", "C"), n_ind, replace = TRUE, prob = runif(3)),
    sites           = sample(paste0("site", 1:4), n_ind, replace = TRUE, prob = runif(4)),
    minimum_n       = sample(2:6, 1),
    minimum_loci    = sample(0:n_loci, 1),
    maf             = runif(1, 0, 0.25),
    max_missingness = runif(1, 0.05, 0.8)
  )
}

n_draws <- 40

test_that("faststats output satisfies invariants for random inputs", {
  set.seed(2024)
  for (draw in seq_len(n_draws)) {
    x <- random_case()
    result <- suppressMessages(faststats(x$gt, x$groups, x$sites,
                                         minimum_n = x$minimum_n,
                                         minimum_loci = x$minimum_loci,
                                         maf = x$maf,
                                         max_missingness = x$max_missingness))
    if (nrow(result) == 0) next

    expect_true(all(c("group", "site", "Ar", "Ho", "He", "uHe", "Fis", "Gst", "Fst",
                      "uFis", "loci", "n") %in% names(result)))
    stat_cols <- setdiff(names(result), c("group", "site"))
    expect_true(all(vapply(result[stat_cols], is.numeric, logical(1))))

    # filters are respected
    expect_true(all(result$n >= x$minimum_n))
    expect_true(all(result$loci >= x$minimum_loci & result$loci <= ncol(x$gt)))
    for (i in seq_len(nrow(result))) {
      in_site <- x$groups == result$group[i] & x$sites == result$site[i]
      expect_equal(result$n[i], sum(in_site))
    }

    # statistics lie in their valid ranges (allowing for rounding to 4 decimal places)
    expect_true(all(result$Ho >= 0 & result$Ho <= 1, na.rm = TRUE))
    expect_true(all(result$He >= 0 & result$He <= 0.5, na.rm = TRUE))
    expect_true(all(result$uHe >= result$He - 1e-4, na.rm = TRUE))
    expect_true(all(result$Ar >= 1 - 1e-4 & result$Ar <= 2 + 1e-4, na.rm = TRUE))
    expect_true(all(result$Fis <= 1 & result$uFis <= 1, na.rm = TRUE))
    expect_true(all(result$Gst <= 1 & result$Fst <= 1, na.rm = TRUE))
  }
})

test_that("per-locus and per-individual statistics satisfy invariants for random inputs", {
  set.seed(7)
  for (draw in seq_len(n_draws)) {
    x <- random_case()
    gt <- x$gt

    hes <- calculate_Hes(gt)
    expect_length(hes, ncol(gt))
    expect_true(all(hes >= 0 & hes <= 0.5, na.rm = TRUE))

    mafs <- get_minor_allele_frequencies(gt)
    expect_true(all(mafs >= 0 & mafs <= 0.5, na.rm = TRUE))
    expect_true(all(get_minor_allele_counts(gt) >= 0))

    ho <- calculate_Ho(gt)
    expect_true(is.na(ho) || (ho >= 0 && ho <= 1))

    ind_ho <- individual_Ho(gt, minimum_loci = 0, maf = 0, max_missingness = 1)
    expect_length(ind_ho, nrow(gt))
    expect_true(all(ind_ho >= 0 & ind_ho <= 1, na.rm = TRUE))

    ar <- calculate_rarified_Ar(gt, x$sites)
    expect_equal(dim(ar), c(length(unique(x$sites)), ncol(gt)))
    expect_true(all(ar >= 1 - 1e-12 & ar <= 2 + 1e-12, na.rm = TRUE))
  }
})

test_that("calculate_Fst is at most 1 and agrees with hierfstat for random inputs", {
  skip_if_not_installed("hierfstat")
  set.seed(11)
  for (draw in seq_len(n_draws)) {
    x <- random_case()
    if (length(unique(x$sites)) < 2) next
    fst <- calculate_Fst(x$gt, x$sites)
    expect_true(all(fst <= 1, na.rm = TRUE))
    expected <- hierfstat::fs.dosage(x$gt, x$sites)$Fs["Fst", names(fst)]
    expect_equal(fst, expected, tolerance = 1e-8)
  }
})

test_that("rarefying a larger sample with the same allele frequencies never gives more alleles", {
  # tripling every individual keeps allele frequencies but triples the sample size; rarefied
  # back to the smaller sample's size, it can only miss alleles the smaller sample has
  set.seed(3)
  for (draw in seq_len(n_draws)) {
    gt <- random_gt(sample(3:20, 1), sample(2:50, 1), missing_rate = 0)
    doubled <- rbind(gt, gt, gt)
    sites <- c(rep("small", nrow(gt)), rep("large", 2 * nrow(gt)))
    ar <- calculate_rarified_Ar(doubled, sites)
    expect_true(all(ar["large", ] <= ar["small", ] + 1e-12))
  }
})

test_that("pi_theta_d_stats output satisfies invariants for random inputs", {
  set.seed(5)
  for (draw in seq_len(n_draws)) {
    x <- random_case()
    result <- suppressWarnings(suppressMessages(
      pi_theta_d_stats(x$gt, x$groups, x$sites, minimum_n = x$minimum_n,
                       minimum_loci = x$minimum_loci, maf = x$maf,
                       max_missingness = x$max_missingness)
    ))
    if (is.null(result)) next

    stat_cols <- setdiff(names(result), c("group", "site"))
    expect_true(all(vapply(result[stat_cols], is.numeric, logical(1))))
    expect_true(all(result$loci >= x$minimum_loci))
    expect_true(all(result$pi >= 0 & result$theta >= 0))
    expect_true(all(result$S <= result$loci))
  }
})
