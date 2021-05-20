
eff_particle_no <- function(w) {
  ess <- 1/sum(w^2)
  return(ess)
}

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


# bootstrap_filter <- function(prior, dyn_model, lik N) {
#   # sample N times from the prior
#   x0 <- prior.sample(N)
#   
#   # importance sampling, zt = \tilde{x_t}
#   for (i in 1:N) {
#     z0[i] <- dyn_model.sample()
#     # update \tilde{x_0:t}
#     _w0[i] <- lik(y0,x0)
#   }
#   # normalise weights
#   w0 <- _w0/sum(_w0)
#   
#   for (t in 1:tt) {
#     
#   }
# }
bootstrap_filter <- function(fk_model, N, tmax, essmin_fn = function(N) N/2) {
  # compute threshold
  essmin <- essmin_fn(N)
  # initialise simulated values of X
  x <- matrix(NA, ncol = N, nrow = (tmax + 1))
  # sample N times from the prior
  x[1, ] <- fk_model$sample_m0(N)
  
  # initialise mean output
  #mx <- rep(NA, tmax)
  
  # initialise weights
  w <- matrix(NA, ncol = N, nrow = (tmax + 1)) # w_t
  W <- matrix(NA, ncol = N, nrow = (tmax + 1)) # W_t
  hw <- matrix(NA, ncol = N, nrow = tmax) # hat{w}_t
  w[1, ] <- exp(fk_model$logG(t, x[1, ]))
  W[1, ] <- w[1, ] / sum(w[1, ])
  
  # initialise ancestor variables
  A <- matrix(NA, ncol = N, nrow = tmax)
  # resampling times
  r <- c()
  
  for (t in 2:(tmax+1)) {
    # resampling step
    if (eff_particle_no(W[t-1, ]) < essmin) {
      r <- c(r, t-1)
      A[t-1, ] <- systematic_resampling(W[t-1, ])
      hw[t-1, ] <- rep(1, N)
    } else {
      A[t-1, ] <- 1:N
      hw[t-1, ] <- w[t-1, ]
    }
    s <- A[t-1, ]
    
    # draw X_t from transition kernel
    x[t, ] <- fk_model$sample_m(x[t-1, s])
    # update mean output
    mx[t-1] <- sum(w * x[t, ])
    # update weights
    w[t, ] <- hw[t-1, ] * exp(fk_model$logG(t, x[t, s]))
    W[t, ] <- w[t, ] / sum(w[t, ])
  }
  return(list(A = A, x = x, hw = hw, w = w, W = W, mx = mx, r = r))
}


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


