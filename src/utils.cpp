#include <Rcpp.h>
#include "utils.h"

Rcpp::NumericVector ncindex(Rcpp::NumericMatrix m, int r, Rcpp::IntegerVector c) {
  int len = c.length();
  Rcpp::NumericVector non_contiguous_row(len);
  
  for (int i=0; i<len; i++) {
    non_contiguous_row[i] = m(r, c(i));
  }
  return non_contiguous_row;
}