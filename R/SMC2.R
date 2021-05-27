
Rcpp::loadModule("particles", TRUE)

# Computation of log(L_T^N) for r_PMMH
log_lik <- function(fk_model, x, A, N, tn) {
  if (is.null(A)) A <- matrix(rep(1:N, tn), nrow = tn, ncol = N, byrow = TRUE)
  if (tn > 1) {
    p1 <- log(sum(exp(fk_model$logG(1, x[1, ]))) / N)
    p2 <- sapply(2:tn, 
                 function(t) log(sum(exp(fk_model$logG(t, x[t, A[t-1, ]]))) / N))
    loglik <- p1 + sum(p2)
  } else {
    loglik <- log(sum(exp(fk_model$logG(1, x[1, ]))) / N)
  }
  return(loglik)
}

# One step particle marginal Metropolis-Hastings
pmmh_onestep_rcpp <- function(fk_model, theta, x, A, mu_prior, sd_prior, sd_prop, ...) {
  # get N and T
  dims <- dim(x)
  N <- dims[2]
  tn <- dims[1] - 1
  # update theta with random walk proposal
  # theta_prop <- theta + rnorm(1, mean = 0, sd = sd_prop)
  theta_prop <- c()
  theta_prop[1] <- theta[1] + rnorm(1, mean = 0, sd = sd_prop[1])
  theta_prop[2] <- exp(log(theta[2]) + rnorm(1, mean = 0, sd = sd_prop[2]))
  # update X and A with boostrap filter
  if (tn > 0) bs_result <- bootstrap_filter_rcpp(fk_model, N, tn)
  else bs_result <- bootstrap_onestep_rcpp(fk_model, N)
  x_prop <- bs_result$x
  A_prop <- bs_result$A
  # compute acceptance probability
  prior_diff <- sum(((theta - mu_prior)^2 - (theta_prop - mu_prior)^2) / (2 * sd_prior^2))
  lr_pmmh <- log_lik(fk_model, x_prop, A_prop, N, tn) - 
    log_lik(fk_model, x, A, N, tn) + prior_diff
  # accept
  if (log(runif(1)) < lr_pmmh) return(list(theta = theta_prop, x = x_prop, A = A_prop))
  # reject
  else return(list(theta = theta, x = x, A = A))
}

# Sequential Monte Carlo squared algorithm
smc_squared_rcpp <- function(Yt, Nx, Nt, sigma, rho, mu_prior, sd_prior, sd_prop, 
                             essmin_fn = function(N) N/2, ...) {
  tmax <- length(Yt) - 1
  # compute threshold
  essmin <- essmin_fn(Nt)
  # initialise thetas
  thetas <- array(NA, dim=c(tmax+1, Nt, 2))
  thetas[1, ,1] <- rnorm(n = Nt, mean = mu_prior[1], sd = sd_prior[1])
  thetas[1, ,2] <- exp(rnorm(n = Nt, mean = log(mu_prior[2]), sd = sd_prior[2]))
  # initialise Nt FK models
  sv_models <- lapply(1:Nt, function(s) new(Bootstrap_SV_C, Yt, thetas[1, s, 1], thetas[1, s, 2], rho))
  # initialise x
  xs <- array(NA, dim = c(tmax+1, Nt, Nx))
  xs[1, , ] <- matrix(unlist(lapply(sv_models, function(x) x$sample_m0(N = Nx))),
                      ncol = Nx, byrow = TRUE)
  # initialise A
  As <- array(NA, dim = c(tmax, Nt, Nx))
  As[1, , ] <- matrix(rep(1:Nx, Nt), nrow = Nt, ncol = Nx, byrow = TRUE)
  # initialise weights
  wm <- matrix(NA, nrow = tmax+1, ncol = Nt) # w^m
  Wm <- matrix(NA, nrow = tmax+1, ncol = Nt) # W^m
  wm[1, ] <- unlist(lapply(
    1:Nt,function(s) log(sum(exp(sv_models[[s]]$logG(1, xs[1, s, ]))) / Nx)))
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
        
        # use population of thetas to improve sampling
        mut <- sapply(1:2, function(k) sum(Wm[t-1,] * thetas[t-1, ,k])/sum(Wm[t-1,]))
        sd_t <- sapply(1:2, function(k) sqrt(sum(Wm[t-1,] * (thetas[t-1,,k]-mut[k])^2)/sum(Wm[t-1,])))
        
        pmmh_results <- pmmh_onestep_rcpp(sv_models[[s]], thetas[t-1, s,], x, A,
                                          mut, sd_prior, sd_t)
        thetas[t, s,] <- pmmh_results$theta
        xs[1:(t-1), s, ] <- pmmh_results$x
        if (!is.null(pmmh_results$A)) As[1:(t-2), s, ] <- pmmh_results$A
      }
      wm[t-1, ] <- 0 # log(1)
      # update Nt FK models
      sv_models <- lapply(1:Nt, 
                          function(s) new(Bootstrap_SV_C, Yt, thetas[t, s, 1], thetas[t, s, 2], rho))
    } else {
      thetas[t, ,] <- thetas[t-1, ,]
    }
    # update ancestor variables
    if (t > 2) As[t-1, , ] <- As[t-2, , ]
    print(colMeans(thetas[t,,]))
    # perform step t for each particle filter
    for (s in 1:Nt) {
      xs[t, s, ] <- sv_models[[s]]$sample_m(xs[t-1, s, As[t-1, s, ]])
      wm[t, s] <- wm[t-1, s] + 
        log(sum(exp(sv_models[[s]]$logG(t, xs[t, s, As[t-1, s, ]]))) / Nx)
    }
    Wm[t, ] <- exp(wm[t, ]) / sum(exp(wm[t, ]))
  }
  return(list(thetas = thetas, xs = xs, As = As, Wm = Wm, ess = ess))
}
