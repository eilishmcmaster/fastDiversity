# Variables used via data.table's `..` prefix and the magrittr `.` placeholder.
# Declared here so R CMD check does not report them as undefined globals.
utils::globalVariables(c(
  ".",
  "..not_missing_loci", "..not_na_loci", "..passing_ac_loci", "..passing_maf_loci"
))
