#include <Rcpp.h>

#ifndef _PARTICLES_H
#define _PARTICLES_H

double eff_particle_no(
    const Rcpp::NumericVector w 
);
// declare systematic resampling function
Rcpp::IntegerVector systematic_resampling(
    const Rcpp::NumericVector W 
);

// declare C++ replica of R class
class Bootstrap_SV_C
{
public:
  Bootstrap_SV_C(Rcpp::NumericVector data, float mu, float sigma, float rho);

  Rcpp::NumericVector sample_m0(int N);
  Rcpp::NumericVector sample_m(Rcpp::NumericVector xp);
  Rcpp::NumericVector logG(int t, Rcpp::NumericVector x);

// private:
  Rcpp::NumericVector _data;
  int tmax;
  float _mu;
  float _sigma;
  float _rho;
  float sigma0;

};


float essmin_fn(int N);

// declare filtering function
Rcpp::List bootstrap_filter(
    Bootstrap_SV_C fk_model, int N, int tmax//, float(*f)(int)
);

Rcpp::List run_bootstrap_filter(Rcpp::NumericVector data, float mu, float sigma, float rho, int N, int tmax);

// declare conversion C++ -> R type



#endif