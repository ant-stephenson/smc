# Computation of L_T^N for r_PMMH
full_lik <- function(fk_model, x, A, theta, N, tn) {
  if (is.null(A)) A <- matrix(rep(1:N, tn), nrow = tn, ncol = N, byrow = TRUE)
  if (tn > 1) {
    p1 <- sum(exp(fk_model$logG(1, x[1, ], theta))) / N
    p2 <- sapply(2:tn, 
                 function(t) sum(exp(fk_model$logG(t, x[t, A[t-1, ]], theta))) / N)
    lik <- p1 * prod(p2)
  } else {
    lik <- sum(exp(fk_model$logG(1, x[1, ], theta))) / N
  }
  return(lik)
}

# Computation of log(L_T^N) for r_PMMH
log_lik <- function(fk_model, x, A, theta, N, tn) {
  if (is.null(A)) A <- matrix(rep(1:N, tn), nrow = tn, ncol = N, byrow = TRUE)
  if (tn > 1) {
    p1 <- log(sum(exp(fk_model$logG(1, x[1, ], theta))) / N)
    p2 <- sapply(2:tn, 
                 function(t) log(sum(exp(fk_model$logG(t, x[t, A[t-1, ]], theta))) / N))
    loglik <- p1 + sum(p2)
  } else {
    loglik <- log(sum(exp(fk_model$logG(1, x[1, ], theta))) / N)
  }
  return(loglik)
}

# One step particle marginal Metropolis-Hastings
pmmh_onestep <- function(fk_model, theta, x, A, mu_prior, sd_prior, sd_prop, ...) {
  # get N and T
  dims <- dim(x)
  N <- dims[2]
  tn <- dims[1] - 1
  # update theta with random walk proposal
  theta_prop <- theta + rnorm(1, mean = 0, sd = sd_prop)
  # update X and A with boostrap filter
  if (tn > 0) bs_result <- bootstrap_filter(fk_model, N, tmax = tn, theta = 0, ...)
  else bs_result <- bootstrap_onestep(fk_model, N, theta = 0, ...)
  x_prop <- bs_result$x
  A_prop <- bs_result$A
  # compute acceptance probability
  #r_num <- dnorm(x = theta_prop, mean = 0, sd = sd_prior, log = TRUE) *
  #  full_lik(fk_model, x_prop, A_prop, theta_prop, N, tn)
  #r_denom <- exp(dnorm(x = theta, mean = 0, sd = sd_prior, log = TRUE)) *
  #  full_lik(fk_model, x, A, theta, N, tn)
  #r_pmmh <- r_num / r_denom
  prior_diff <- ((theta - mu_prior)^2 - (theta_prop - mu_prior)^2) / (2 * sd_prior^2)
  lr_pmmh <- #dnorm(x = theta_prop, mean = mu_prior, sd = sd_prior, log = TRUE) +
    log_lik(fk_model, x_prop, A_prop, theta_prop, N, tn) -
  #  dnorm(x = theta_prop, mean = mu_prior, sd = sd_prior, log = TRUE) -
    log_lik(fk_model, x, A, theta, N, tn) + prior_diff
  # accept
  if (log(runif(1)) < lr_pmmh) return(list(theta = theta_prop, x = x_prop, A = A_prop))
  # reject
  else return(list(theta = theta, x = x, A = A))
}

# Sequential Monte Carlo squared algorithm
smc_squared <- function(Yt, Nx, Nt, sigma, rho, mu_prior, sd_prior, sd_prop, 
                        essmin_fn = function(N) N/2, ...) {
  tmax <- length(Yt) - 1
  # compute threshold
  essmin <- essmin_fn(Nt)
  # initialise thetas
  thetas <- matrix(NA, nrow = (tmax + 1), ncol = Nt)
  thetas[1, ] <- rnorm(n = Nt, mean = mu_prior, sd = sd_prior)
  # initialise Nt FK models
  sv_models <- lapply(1:Nt,
                      function(s) Bootstrap_SV$new(data = Yt, mu = thetas[1, s],
                                                   sigma = sigma, rho = rho))
  # initialise x
  xs <- array(NA, dim = c(tmax+1, Nt, Nx))
  xs[1, , ] <- matrix(unlist(lapply(sv_models,
                                    function(x) x$sample_m0(N = Nx))),
                      ncol = Nx, byrow = TRUE)
  # initialise A
  As <- array(NA, dim = c(tmax, Nt, Nx))
  As[1, , ] <- matrix(rep(1:Nx, Nt), nrow = Nt, ncol = Nx, byrow = TRUE)
  # initialise weights
  wm <- matrix(NA, nrow = tmax+1, ncol = Nt) # w^m
  Wm <- matrix(NA, nrow = tmax+1, ncol = Nt) # W^m
  wm[1, ] <- unlist(lapply(1:Nt,
                           function(s) {
                             log(sum(exp(sv_models[[s]]$logG(1, xs[1, s, ], 
                                                             thetas[1, s]))) / Nx)
                           }))
  Wm[1, ] <- exp(wm[1, ]) / sum(exp(wm[1, ]))
  # initialise ESS vector
  ess <- c()
  # SMC^2 algorithm
  for (t in 2:(tmax+1)) {
    # compute ess
    ess[t-1] <- eff_particle_no(Wm[t-1, ])
    print(ess[t-1])
    if (ess[t-1] < essmin) {
      # move particles through PMMH kernel
      for (s in 1:Nt) {
        if (t == 2) A <- NULL
        else A <- matrix(As[1:(t-2), s, ], nrow = t-2)
        x <- matrix(xs[1:(t-1), s, ], nrow = t-1)
        pmmh_results <- pmmh_onestep(sv_models[[s]], thetas[t-1, s], x, A, 
                                     mu_prior, sd_prior, sd_prop, ...)
        thetas[t, s] <- pmmh_results$theta
        xs[1:(t-1), s, ] <- pmmh_results$x
        if (!is.null(pmmh_results$A)) As[1:(t-2), s, ] <- pmmh_results$A
      }
      wm[t-1, ] <- 0 # log(1)
      # update Nt FK models
      sv_models <- lapply(1:Nt, 
                          function(s) Bootstrap_SV$new(data = Yt, mu = thetas[t, s],
                                                       sigma = sigma, rho = rho))
    } else {
      thetas[t, ] <- thetas[t-1, ]
    }
    # update ancestor variables
    if (t > 2) As[t-1, , ] <- As[t-2, , ]
    print(mean(thetas[t, ]))
    # perform step t for each particle filter
    xs[t, , ] <- matrix(unlist(lapply(1:Nt, function(s) {
      sv_models[[s]]$sample_m(xs[t-1, s, As[t-1, s, ]])
    })), ncol = Nx, byrow = TRUE)
    wm[t, ] <- wm[t-1, ] + unlist(lapply(1:Nt,
                                         function(s) {
                                           log(sum(exp(sv_models[[s]]$logG(
                                             t, xs[t, s, As[t-1, s, ]],
                                             thetas[t, s]))) / Nx)
                                         }))
    Wm[t, ] <- exp(wm[t, ]) / sum(exp(wm[t, ]))
  }
  return(list(thetas = thetas, xs = xs, As = As, wm = wm, Wm = Wm, ess = ess))
}
