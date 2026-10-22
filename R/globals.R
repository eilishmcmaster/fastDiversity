# Variables used via data.table's `..` prefix, the magrittr `.` placeholder,
# and functions created at runtime by Rcpp::cppFunction(). Declared here so
# R CMD check does not report them as undefined globals.
utils::globalVariables(c(
  ".",
  "..not_missing_loci", "..not_na_loci", "..passing_ac_loci", "..passing_maf_loci",
  "calculate_Ho_cpp", "calculate_He_cpp",
  "calculate_Ho_vec_cpp", "calculate_Hes_vec_cpp"
))
