#include <Rcpp.h>
#include "utils.h"
#include "particles.h"
using namespace Rcpp;

double eff_particle_no(
    const NumericVector w 
) 
{
  double ess = 1.0/sum(w * w);
  return ess;
}

// [[Rcpp::export(name = "systematic_resampling_rcpp")]]
IntegerVector systematic_resampling(const NumericVector W) {
  int N = W.length();
  NumericVector v = cumsum(W);
  float s = R::runif(0, 1) * 1.0/(float)N;
  IntegerVector A(N);
  int m = 0;
  for (int n=0; n<N; n++) {
    while(v[n] > s) {
      m += 1;
      s += 1.0/(float)N;
    }
    A[n] = (m == 0) ? 1 : m;
  }
  return A;
}

Bootstrap_SV_C::Bootstrap_SV_C(NumericVector data, float mu, float sigma, float rho)
{
  this->_data = data;
  this->_mu = mu;
  this->_sigma = sigma;
  this->_rho = rho;
  this->tmax = data.length();
  this->sigma0 = sigma / sqrt(1.0 - rho * rho);
}

NumericVector Bootstrap_SV_C::sample_m0(int N)
{
  return rnorm(N, this->_mu, this->sigma0);
}

NumericVector Bootstrap_SV_C::sample_m(NumericVector xp)
{
  int N = xp.length();
  NumericVector sample(N);

  for (int i=0; i<N; i++) {
    sample[i] = R::rnorm(this->_mu + this->_rho * (xp[i] - this->_mu), this->_sigma);
  }
  return sample;
}

NumericVector Bootstrap_SV_C::logG(int t, NumericVector x)
{
  int N = x.length();
  NumericVector logg(N);

  for (int i=0; i<N; i++) {
    logg[i] = R::dnorm(this->_data[t], 0.0, sqrt(exp(x[i])), true);
  }
  return logg;
}


float essmin_fn(int N) {
  return (float)N/2.0;
}


List bootstrap_filter_rcpp(Bootstrap_SV_C fk_model, int N, int tmax) {//, float(*f)(int) = [](int N) {return essmin_fn(N);}) {
  //float essmin = (*f)(N);
  float essmin = essmin_fn(N);
  
  // initialise simulated values of X
  NumericMatrix x(tmax+1, N);
  // sample N times from the prior
  x(0, _) = fk_model.sample_m0(N);
  
  // initialise mean and sd output
  NumericVector mx(tmax+1);
  NumericVector sdx(tmax+1);
  
  // initialise ess storage
  IntegerVector ess(tmax);
  
  // initialise weights
  NumericVector w(N);
  NumericMatrix W(tmax+1, N);
  NumericVector hw(N);
  
  w = exp(fk_model.logG(1, x(0, _)));
  W(0, _) = w / sum(w);
  
  // initialise ancestor variables
  IntegerMatrix A(tmax, N);
  
  // resampling times
  IntegerVector r;
  
  for (int t=1; t<=tmax; t++) {
    ess(t-1) = eff_particle_no(W(t-1, _));
    if (ess(t-1) < essmin) {
      r.push_back(t-1);
      A(t-1, _) = systematic_resampling(W(t-1, _));
      hw = 0*hw + 1;
    }
    else {
      for (int i=0; i<N; i++) {
        A(t-1, i) = i + 1;
      }
      hw = w;
    }
    IntegerMatrix::Row s = A(t-1, _);
    s = s - 1;
    
    // draw X_t from transition kernel
    NumericVector xp = ncindex(x, t-1, s);
    x(t, _) = fk_model.sample_m(xp);
    
    // update weights
    NumericVector xt = ncindex(x, t, s);
    w = hw * exp(fk_model.logG(t, xt));
    W(t, _) = w / sum(w);
    
    // update mean and sd output
    mx(t) = sum(W(t, _) * x(t, _));
    sdx(t) = sqrt(sum(pow(x(t, _) - mx(t), 2.0) / (float)(N-1)));
  }
  
  List output = List::create(_["A"] = A, _["x"] = x, _["hw"] = hw, _["W"] = W, 
                             _["mx"] = mx, _["sdx"] = sdx, _["r"] = r, _["ess"] = ess);
  
  return output;
}


RCPP_EXPOSED_CLASS(Bootstrap_SV_C);
RCPP_MODULE(particles) {
  
  Rcpp::class_<Bootstrap_SV_C>("Bootstrap_SV_C")
  
  .constructor<Rcpp::NumericVector, float, float, float>()
  
  .field("_data", &Bootstrap_SV_C::_data)
  .field("tmax", &Bootstrap_SV_C::tmax)
  .field("_mu", &Bootstrap_SV_C::_mu)
  .field("_sigma", &Bootstrap_SV_C::_sigma)
  .field("_rho", &Bootstrap_SV_C::_rho)
  .field("sigma0", &Bootstrap_SV_C::sigma0)
  
  .method("sample_m0", &Bootstrap_SV_C::sample_m0)
  .method("sample_m", &Bootstrap_SV_C::sample_m)
  .method("logG", &Bootstrap_SV_C::logG)
  ;
  
  function("bootstrap_filter_rcpp", &bootstrap_filter_rcpp);
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
# set.seed(1)
# 
# tmax <- 100
# mu <- -1
# rho <- 0.95
# sigma <- 0.15
# 
# N <- 1000
# 
# Xt <- generate_SV_data(mu, rho, sigma, tmax)
# Yt <- as.matrix(rnorm(tmax, mean = 0, sd = sqrt(exp(Xt))))
# 
# boot_sv <- new(Bootstrap_SV_C, Yt, mu, sigma, rho)
# 
# output <- bootstrap_filter_rcpp(boot_sv, N, tmax)
*/
