
eff_particle_no <- function(w) {
  ess <- 1/sum(w^2)
  return(ess)
}

systematic_resampling <- function(w) {
  U <- runif(1)
  N <- length(w)
  v <- cumsum(N * w)
  s <- U
  m <- 1
  A <- 0*1:N
  for (n in 1:N) {
    while (v[m] < s) {
      m <- m + 1
      A[n] <- m
      s <- s + 1
    }
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
bootstrap_filter <- function(fk_model, N, tt, ESSmin=function(N) N/2) {
  ## initialise
  # sample N times from the prior
  x <- matrix(rep(NA, (tt+1)*N), nrow=(tt+1))
  x[1, ] <- fk_model$sample_m0(N)
  
  # mean output
  mx <- rep(NA, tt)
  
  # initialise weights
  w <- rep(1/N, N)
  s <- 1:N

  ## run for t:T
  for (t in 2:(tt+1)) {
    # print(eff_particle_no(w))
    # print(ESSmin(N))
    # if (any(is.na(w))) {print(w)}
    # browser()
    
    
    x[t, ] <- fk_model$sample_m(t, x[t-1, ])
    wt <- exp(fk_model$logG(t, x[t, s], x[t, ]))
    
    w <- wt/sum(wt)
    
    mx[t-1] <- sum(w*x[t,])
    
    if (eff_particle_no(w) < ESSmin(N)) {
      s <- systematic_resampling(w)
      w <- rep(1/N, N)
    }

    s <- sample(1:N, size=N, replace=TRUE, prob=w)
    
    x <- x[, s]
  }
  return(list(w=w, x=x, mx=mx))
}


require(methods)
Bootstrap_SV <- setRefClass("Bootstrap_SV",
                            fields=list(data="matrix", mu="numeric", sigma="numeric", rho="numeric", tt="numeric", sigma0="numeric"),
                            methods =list(
                              initialize = function(data, mu=0, sigma=1, rho=0.95) {
                                .self$data = data
                                .self$tt = length(data)
                                .self$mu = mu
                                .self$sigma = sigma
                                .self$rho = rho
                                .self$sigma0 = sigma/sqrt(1-rho^2)
                              },
                              sample_m0 = function(N) {
                                rnorm(N, mean=.self$mu, sd=.self$sigma0)
                              },
                              sample_m = function(t, xp) {
                                rnorm(length(xp), mean=.self$mu + .self$rho * (xp - .self$mu), sd = .self$sigma)
                              },
                              logG = function(t, xp, x) {
                                dnorm(.self$data[t], mean=xp, sd=sqrt(exp(0.5 * x)), log=TRUE)
                              }
                            ))


