#include <Rcpp.h>

#ifndef _PARTICLES_H
#define _PARTICLES_H

double eff_particle_no(
    const Rcpp::NumericVector w 
);
Rcpp::IntegerVector systematic_resampling(
    const Rcpp::NumericVector W 
);

class Bootstrap_SV
{
public:
  Bootstrap_SV(Rcpp::NumericVector data, float mu, float sigma, float rho);
  
  Rcpp::NumericVector sample_m0(int N);
  Rcpp::NumericVector sample_m(Rcpp::NumericVector xp);
  Rcpp::NumericVector logG(int t, Rcpp::NumericVector x);
  
private:
  Rcpp::NumericVector _data;
  int tmax;
  float _mu;
  float _sigma;
  float _rho;
  float sigma0;
  
};

struct FilterOutput
{
  Rcpp::IntegerMatrix A;
  Rcpp::NumericMatrix x;
  Rcpp::NumericVector hw;
  Rcpp::NumericMatrix W;
  Rcpp::NumericVector mx;
  Rcpp::IntegerVector r;
  
  FilterOutput(const Rcpp::IntegerMatrix A, const Rcpp::NumericMatrix x, 
               const Rcpp::NumericVector hw, const Rcpp::NumericMatrix W,
               const Rcpp::NumericVector mx, const Rcpp::IntegerVector r) 
    : A(A), x(x), hw(hw), W(W), mx(mx), r(r) {}
};

float essmin_fn(int N);

FilterOutput bootstrap_filter(
    Bootstrap_SV fk_model, int N, int tmax, float(*f)(int)
);

#endif