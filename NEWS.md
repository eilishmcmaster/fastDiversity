# fastDiversity (development version)

## Changes to statistics
- Documentation now explains how monomorphic sites affect each statistic,
  including that `pi` and `theta` from `pi_theta_d_stats()` are sums over
  segregating sites rather than per-site values.
- `faststats()`: the previous `Fst` column, (Ht - Hs) / Ht, was mislabelled as
  Weir & Cockerham (1984). It is now named `Gst` (Nei 1973). A new `Fst`
  column gives population-specific FST (Weir & Goudet 2017) for each site,
  relative to the other retained sites in the genetic group.
- New `calculate_Fst()`: population-specific FST (Weir & Goudet 2017); matches
  `hierfstat::fs.dosage()`.
- `faststats()`: `Ar` is now rarefied allelic richness, as documented, using
  `calculate_rarified_Ar()` across the retained sites in each genetic group.
  Previously it was the mean number of alleles per individual (1 + Ho).
- `calculate_Ar()`: now returns the (unrarefied) mean number of alleles
  observed per locus. Previously it returned 1 + Ho.
- `calculate_rarified_Ar()`: uses log-scale binomial coefficients (previously
  overflowed to 0 for large samples), returns NA for loci not genotyped at
  every site, returns site and locus names, and handles single-individual
  sites; matches `hierfstat::allelic.richness()` on complete data.

## Bug fixes
- `make_allele_list()` / `matcher()`: alleles were assigned using DArT 1-row
  coding (2 = heterozygote) instead of the package's 0=aa, 1=aA, 2=AA coding.
  Private and total allele counts from `calculate_private_alleles()` change
  as a result, and now agree with `calculate_total_alleles()`.
- `faststats()`: statistic columns are now numeric (previously character).
- `faststats()`: `return_locus_stats = TRUE` no longer errors when
  `run_HWE_test = FALSE`; it now runs the HWE test automatically.
- `resampling_loci_CI()`, `resampling_individuals_CI()`: no longer call
  `set.seed(14)` internally, which overwrote the user's random number state.
  Use `set.seed()` before calling them for reproducible intervals.
- `MC_HWEtest_of_loci()`: matrix input now works (previously only data
  frames), and the per-locus "genotypes are not labelled" warnings are gone.
- `individual_Ho()`: no longer errors when a group has too few loci.
- `faststats()`: statistics are calculated at full precision and rounded only
  in the output (per-locus He was previously rounded before averaging).
- `faststats()`: groups where no loci pass the missingness filter are skipped
  with a message instead of erroring.
- `faststats()`: bootstrap CIs no longer error when all resampled values are NA.
- `pi_theta_d_stats()`: `minimum_loci` is now applied (it was ignored),
  statistic columns are numeric, and sites with a single remaining locus or no
  segregating sites no longer error.

## Other changes
- C++ helpers are compiled with the package instead of at runtime with
  `Rcpp::cppFunction()`, so bootstrapping is much faster and no compiler is
  needed after installation.
- Progress output uses `message()` instead of `print()`/`cat()`, so it can be
  silenced with `suppressMessages()`.
- Unused `future.apply` dependency removed.
- Example dataset documented (`?example_data`).
- Expanded test suite with checks against hand-calculated values.
- pkgdown documentation site and R CMD check GitHub Action added.

# fastDiversity 1.3.0

## Changes to existing functions
- `calculate_pi_theta_d()`: reimplemented with a vectorised approach for 
  improved performance. Non-segregating and monomorphic sites are now removed 
  prior to calculation. The function now accepts the same `gt` matrix format 
  used throughout the package. Return value structure has changed: the function 
  now returns `tajima_d`, `num_sites`, `raw_pi`, `watterson_theta`, and 
  `d_stdev` (previously returned `pi`, `S`, `theta`, `D`, and `per_site_stats`).