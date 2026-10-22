#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector calculate_Ho_vec_cpp(NumericMatrix gt_mat) {
  int n_col = gt_mat.ncol();
  int n_row = gt_mat.nrow();
  NumericVector ho(n_col);
  
  for (int j = 0; j < n_col; ++j) {
    int heteros = 0;
    int total = 0;
    for (int i = 0; i < n_row; ++i) {
      if (!NumericVector::is_na(gt_mat(i,j))) {
        total++;
        if (gt_mat(i,j) == 1) heteros++;
      }
    }
    ho[j] = total > 0 ? (double)heteros / total : NA_REAL;
  }
  return ho;
}

// [[Rcpp::export]]
NumericVector calculate_Hes_vec_cpp(NumericMatrix gt_mat) {
  int n_col = gt_mat.ncol();
  int n_row = gt_mat.nrow();
  NumericVector he(n_col);
  
  for (int j = 0; j < n_col; ++j) {
    int allele_count = 0;
    int alt_count = 0;
    for (int i = 0; i < n_row; ++i) {
      if (!NumericVector::is_na(gt_mat(i,j))) {
        allele_count += 2;
        alt_count += gt_mat(i,j);
      }
    }
    if (allele_count > 0) {
      double alt_freq = (double)alt_count / allele_count;
      he[j] = 2 * alt_freq * (1 - alt_freq);
    } else {
      he[j] = NA_REAL;
    }
  }
  return he;
}
