#include <Rcpp.h>
#include "utils.h"
using namespace Rcpp;

NumericVector index(NumericMatrix m, int r, IntegerVector c) {
  int len = c.length();
  NumericVector non_contiguous_row(len);
  
  for (int i=0; i<=len; i++) {
    non_contiguous_row[i] = m(r, c(i));
  }
  return non_contiguous_row;
}