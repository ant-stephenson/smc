#' Effective Particle Number
#' 
#' Uses the weight vector input to compute an estimate of effective particle number for use in particle filters.
#' @name eff_particle_no
#' @export eff_particle_no
#'
#' @field w vector of weights
# eff_particle_no <- function(w) {
#   ess <- 1/sum(w^2)
#   return(ess)
# }

#' Systematic resampling algorithm
#' 
#' Samples a single uniform random variable U and then assigns values U(n) = (n - 1 + U)/N.
#' Used to generate a vector of integers to index the particles.
#' @name systematic_resampling
#' @export systematic_resampling
systematic_resampling <- function(W) {
  N <- length(W)
  v <- cumsum(W)
  s <- runif(1, 0, 1/N)
  A <- rep(0, N)
  m <- 0
  for (n in 1:N) {
    while (v[n] > s) {
      m <- m + 1
      s <- s + 1/N
    }
    A[n] <- ifelse(m == 0, 1, m)
  }
  return(A)
}

#' Bootstrap filter
#' 
#' Implements a bootstrap particle filter.
#' Takes an object of class Bootstrap_SV_C as argument, with number of particles
#' N and number of time steps to run.
#' @name bootstrap_filter_rcpp
#' @export bootstrap_filter_rcpp
#' @field fk_model Bootstrap_SV_C object
#' @field N number of particles
#' @field tmax number of steps
bootstrap_filter <- function(fk_model, N, tmax,
                            resampling = systematic_resampling,
                            essmin_fn = function(N) N/2) {
  # compute threshold
  essmin <- essmin_fn(N)
  # initialise simulated values of X
  x <- matrix(NA, nrow = (tmax + 1), ncol = N)
  # sample N times from the prior
  x[1, ] <- fk_model$sample_m0(N)
  # initialise vector of ess
  ess <- c()
  
  # initialise weights
  w <- matrix(NA, nrow = (tmax + 1), ncol = N) # w_t
  W <- matrix(NA, nrow = (tmax + 1), ncol = N) # W_t
  hw <- matrix(NA, nrow = tmax, ncol = N) # hat{w}_t
  w[1, ] <- fk_model$logG(1, x[1, ])
  W[1, ] <- exp(w[1, ]) / sum(exp(w[1, ]))
  
  # initialise mean and sd output
  mx <- c()
  sdx <- c()
  mx[1] <- sum(W[1, ] * x[1, ])
  sdx[1] <- sum((x[1, ] - mx[1])^2) / (N-1)
  
  # initialise ancestor variables
  A <- matrix(NA, nrow = tmax, ncol = N)
  # resampling times
  r <- c()
  
  for (t in 2:(tmax+1)) {
    # compute ess
    ess[t-1] <- eff_particle_no(W[t-1, ])
    # resampling step
    if (ess[t-1] < essmin) {
      r <- c(r, t-1)
      A[t-1, ] <- resampling(W[t-1, ])
      hw[t-1, ] <- rep(0, N)
    } else {
      A[t-1, ] <- 1:N
      hw[t-1, ] <- w[t-1, ]
    }
    # draw X_t from transition kernel
    x[t, ] <- fk_model$sample_m(x[t-1, A[t-1, ]])
    # update weights
    w[t, ] <- hw[t-1, ] + fk_model$logG(t, x[t, A[t-1, ]])
    W[t, ] <- exp(w[t, ]) / sum(exp(w[t, ]))
    # update mean and sd output
    mx[t] <- sum(W[t, ] * x[t, ])
    sdx[t] <- sqrt(sum((x[t, ] - mx[t])^2) / (N-1))
  }
  return(list(A = A, x = x, hw = hw, w = w, W = W, 
              mx = mx, sdx = sdx, r = r, ess = ess))
}

bootstrap_onestep <- function(fk_model, N, theta = NULL,
                             resampling = systematic_resampling,
                             essmin_fn = function(N) N/2) {
  # sample N times from the prior of X
  x <- matrix(fk_model$sample_m0(N), ncol = N, nrow = 1)
  # initialise weights
  w <- matrix(fk_model$logG(1, x[1, ]), ncol = N, nrow = 1) # w_t
  W <- matrix(exp(w[1, ]) / sum(exp(w[1, ])), ncol = N, nrow = 1) # W_t
  if (eff_particle_no(W) < essmin_fn(N)) {
    A <- matrix(resampling(W), ncol = N, nrow = 1)
  } else {
    A <- matrix(1:N, ncol = N, nrow = 1)
  }
  # draw X_t from transition kernel
  x <- matrix(fk_model$sample_m(x[, A[1, ]]), ncol = N, nrow = 1)
  return(list(A = A, x = x))
}

#' Bootstrap class for Stochastic Volatility (SV) model
#' 
#' Initialised with a data vector of observations (Yt), an estimate mean and 
#' standard deviation parameter and rho controls the autocorrelation.
#' 
#' The object defines methods to sample from the initial distribution M0 and 
#' the subsequent transition kernels Mt and calculate the log probability from 
#' the potential function Gt.
#' @name Bootstrap_SV_C
#' @export Bootstrap_SV_C
#' @exportClass Bootstrap_SV_C
#' @field data vector of observations
#' @field mu estimate of mean of latent var (float)
#' @field sigma estimate of sd of latent
#' @field rho autocorrelation parameter
require(methods)
Bootstrap_SV <- setRefClass("Bootstrap_SV",
                             fields = list(
                               data = "matrix", 
                               mu = "numeric", 
                               sigma = "numeric", 
                               rho = "numeric", 
                               tmax = "numeric", 
                               sigma0 = "numeric"),
                             methods = list(
                               initialize = function(data, mu = 0, sigma = 1, rho = 0.95) {
                                 .self$data = data
                                 .self$tmax = length(data)
                                 .self$mu = mu
                                 .self$sigma = sigma
                                 .self$rho = rho
                                 .self$sigma0 = sigma / sqrt(1 - rho^2)
                               },
                               sample_m0 = function(N) {
                                 rnorm(N, mean = .self$mu, sd = .self$sigma0)
                               },
                               sample_m = function(xp) {
                                 rnorm(length(xp), 
                                       mean = .self$mu + .self$rho * (xp - .self$mu), 
                                       sd = .self$sigma)
                               },
                               logG = function(t, x) {
                                 dnorm(.self$data[t], mean = 0, 
                                       sd = sqrt(exp(x)), log = TRUE)
                               }
                             ))


