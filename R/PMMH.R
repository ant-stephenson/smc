# Computation of L_T^N for r_PMMH
full_lik <- function(fk_model, x, A, theta) {
  # get N and T
  dims <- dim(A)
  N <- dims[1]
  tmax <- dims[2]
  # compute the likelihood
  p1 <- sum(exp(fk_model$logG(t = 1, x = x[1, ], mean = theta))) / N
  p2 <- sapply(1:tmax, function(t) {
    sum(exp(fk_model$logG(t = A[t, ], x = x[t+1, ], mean = theta))) / N
  })
  lik <- p1 * prod(p2)
  return(lik)
}

# One step particle marginal Metropolis hastings
pmmh_onestep <- function(fk_model, theta, x, A, sd_prop, sd_prior, ...) {
  # update theta with random walk proposal
  theta_prop <- theta + rnorm(1, mean = 0, sd = sd_prop)
  # update X and A with bootstrap filter
  bs_result <- bootstrap_filter(fk_model = fk_model, N = N, tmax = tmax, ...)
  x_prop <- bs_result$x
  A_prop <- bs_result$A
  # remove object to save memory
  rm(list = bs_result, envir = .GlobalEnv)
  # compute acceptance probability
  r_num <- dnorm(x = theta_prop, mean = 0, sd = sd_prior) *
    full_lik(fk_model, x_prop, A, theta_prop)
  r_denom <- dnorm(x = theta, mean = 0, sd = sd_prior) *
    full_lik(fk_model, x, A, theta)
  r_pmmh <- r_num / r_denom
  if(runif(1) < r_pmmh) return(list(theta = theta_prop, x = x_prop, A = A_prop))
  else return(list(theta = theta, x = x, A = A))
}

smc_squared <- function(Yt, Nx, Nt, sigma, rho, sd_prior, sd_prop, 
                        essmin_fn = function(N) N/2, ...) {
  tmax <- length(Yt) - 1
  # compute threshold
  essmin <- essmin_fn(Nt)
  # initialise thetas
  thetas <- matrix(NA, nrow = (tmax + 1), ncol = Nt)
  thetas[1, ] <- rnorm(n = Nt, mean = 0, sd = sd_prior)
  # initialise FK models
  sv_models <- lapply(1:Nt, function(s) Bootstrap_SV$new(data = Yt, mu = thetas[1, s], 
                                                         sigma = sigma, rho = rho))
  # initialise x
  xs <- array(NA, dim = c(tmax+1, Nt, Nx))
  xs[1, , ] <- unlist(lapply(sv_models, function(x) x$sample_m0(N = Nx))) %>% 
    matrix(ncol = Nx) %>% t()
  # initialise A
  As <- array(NA, dim = c(tmax+1, Nt, tmax+1))
  As[1, , ] <- matrix(1:(tmax+1), nrow = Nt, ncol = (tmax + 1), byrow = TRUE)
  # initialise weights
  wm <- matrix(NA, nrow = Nx, ncol = Nt) # w^m
  Wm <- matrix(NA, nrow = Nx, ncol = Nt) # W^m
  wm[1, ] <- unlist(lapply(1:Nt, 
                           function(s) sum(exp(sv_models[[s]]$logG(1, xs[1, s, ]))) / Nx))
  Wm[1, ] <- wm[1, ] / sum(wm[1, ])
  # init ess vector
  ess <- c()
  # SMC^2 steps
  for (t in 2:(tmax+1)) {
    # compute ess
    ess[t-1] <- eff_particle_no(Wm[t-1, ])
    print(ess[t-1])
    if (ess[t-1] < essmin) {
      # move particles through PMMH kernel
      for (s in 1:Nt) {
        pmmh_results <- pmmh_onestep(fk_model = sv_models[[s]], theta = thetas[t-1, ], 
                                     x = xs[t-1, s, ], A = As[t-1, s, ], 
                                     sd_prop = sd_prop, sd_prior = sd_prior, ...)
        thetas[t, ] <- pmmh_results$theta
        xs[t, s, ] <- pmmh_results$x
        As[t, s, ] <- pmmh_results$A
      }
      wm[t-1, ] <- 1
    } else {
      thetas[t, ] <- thetas[t-1, ]
    }
    # update FK models - rearrange Yt according to A_t^m
    sv_models <- lapply(1:Nt, 
                        function(s) Bootstrap_SV$new(data = as.matrix(Yt[As[t, s, ]]), 
                                                     mu = thetas[t, s], 
                                                     sigma = sigma, rho = rho))
    # resample x
    xs[t, , ] <- unlist(lapply(sv_models, 
                               function(x) x$sample_m(xp = xs[t-1, s, ]))) %>% 
      matrix(ncol = Nx) %>% t()
    # update weights
    wm[t, ] <- wm[t-1, ] * 
      unlist(lapply(1:Nt, function(s) sum(exp(sv_models[[s]]$logG(t, xs[t, s, ]))) / Nx))
    Wm[t, ] <- wm[t, ] / sum(wm[t, ])
  }
  return(list(thetas = thetas, xs = xs, As = As, wm = wm, Wm = Wm, ess = ess))
}
