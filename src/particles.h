#include <Rcpp.h>
// #include <RcppCommon.h>

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
class Bootstrap_SV
{
public:
  Bootstrap_SV(Rcpp::NumericVector data, float mu, float sigma, float rho);
  
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
    Bootstrap_SV fk_model, int N, int tmax//, float(*f)(int)
);

Rcpp::List run_bootstrap_filter(Rcpp::NumericVector data, float mu, float sigma, float rho, int N, int tmax);

// declare conversion C++ -> R type
// namespace Rcpp {
//   template <> SEXP wrap(const Bootstrap_SV& obj);
// }

RCPP_MODULE(filter_module) {
  
  Rcpp::class_<Bootstrap_SV>("Bootstrap_SV")
  
  .constructor<Rcpp::NumericVector, float, float, float>()
  
  .field("_data", &Bootstrap_SV::_data)
  .field("tmax", &Bootstrap_SV::tmax)
  .field("_mu", &Bootstrap_SV::_mu)
  .field("_sigma", &Bootstrap_SV::_sigma)
  .field("_rho", &Bootstrap_SV::_rho)
  .field("sigma0", &Bootstrap_SV::sigma0)
  
  .method("sample_m0", &Bootstrap_SV::sample_m0)
  .method("sample_m", &Bootstrap_SV::sample_m)
  .method("logG", &Bootstrap_SV::logG)
  ;
}

#endif