# Small genotype matrix with hand-calculated expected values, used across tests.
#
#        L1  L2  L3
#   i1    0   0   2
#   i2    1   0   2
#   i3    2   1   2
#   i4    1  NA   1
#
# L1: alt alleles 4/8, p = 1/2, He = 2(1/2)(1/2) = 1/2,      Ho = 2/4
# L2: alt alleles 1/6, p = 1/6, He = 2(1/6)(5/6) = 10/36,    Ho = 1/3
# L3: alt alleles 7/8, p = 7/8, He = 2(7/8)(1/8) = 14/64,    Ho = 1/4
toy_gt <- matrix(
  c(0, 1, 2, 1,
    0, 0, 1, NA,
    2, 2, 2, 1),
  nrow = 4,
  dimnames = list(paste0("i", 1:4), paste0("L", 1:3))
)

toy_Hes <- c(1/2, 10/36, 14/64)
toy_ns  <- c(4, 3, 4)
