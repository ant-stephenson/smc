# Computation of L_T^N for r_PMMH
full_lik <- function(fk_model, x, A, theta, N, tn) {
  if(is.null(A)) A <- matrix(rep(1:N, tn), nrow = tn, ncol = N, byrow = TRUE)
  # compute the likelihood
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

# One step particle marginal Metropolis hastings
pmmh_onestep <- function(fk_model, theta, x, A, sd_prop, sd_prior, ...) {
  # get N and T
  dims <- dim(x)
  N <- dims[2]
  tn <- dims[1] - 1
  # update theta with random walk proposal
  theta_prop <- theta + rnorm(1, mean = 0, sd = sd_prop)
  # update X and A with bootstrap filter
  if (tn > 0) bs_result <- bootstrap_filter(fk_model = fk_model, N = N, tmax = tn, ...)
  else bs_result <- bootstrap_onestep(fk_model, N = N, ...)
  x_prop <- bs_result$x
  A_prop <- bs_result$A
  # compute acceptance probability
  r_num <- dnorm(x = theta_prop, mean = 0, sd = sd_prior) *
    full_lik(fk_model, x_prop, A, theta_prop, N, tn)
  r_denom <- dnorm(x = theta, mean = 0, sd = sd_prior) *
    full_lik(fk_model, x, A, theta, N, tn)
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
  xs[1, , ] <- t(matrix(unlist(lapply(sv_models, function(x) x$sample_m0(N = Nx))), 
                        ncol = Nx))
  # initialise A
  As <- array(NA, dim = c(tmax+1, Nt, Nx))
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
                                     x = t(as.matrix(xs[1:(t-1), s, ])), 
                                     A = ifelse(t > 2, t(as.matrix(As[1:(t-2), s, ])), 
                                                NULL), 
                                     sd_prop = sd_prop, sd_prior = sd_prior, ...)
        thetas[t, ] <- pmmh_results$theta
        xs[t, s, ] <- pmmh_results$x
        As[t-1, s, ] <- pmmh_results$A
        # update FK models - rearrange Yt according to A_t^m
        sv_models <- lapply(1:Nt, 
                            function(s) Bootstrap_SV$new(data = Yt, 
                                                         mu = thetas[t, s], 
                                                         sigma = sigma, rho = rho))
      }
      wm[t-1, ] <- 1
    } else {
      thetas[t, ] <- thetas[t-1, ]
      if (t == 2) {
        # make an initial matrix for ancestor variables
        As[t-1, , ] <- matrix(rep(1:Nx, tmax), nrow = tmax, ncol = Nx, byrow = TRUE)
      } else {
        As[t-1, , ] <- As[t-2, , ]
      }
    }
    # resample x
    xs[t, , ] <- t(matrix(unlist(lapply(1:Nt, 
                                        function(s) Wm[t-1, s] * 
                                          sv_models[[s]]$sample_m(xs[t-1, s, As[t-1, s, ]]))), 
                          ncol = Nx))
    # update weights
    wm[t, ] <- wm[t-1, ] * 
      unlist(lapply(1:Nt, function(s) sum(exp(sv_models[[s]]$logG(t, xs[t, s, ]))) / Nx))
    Wm[t, ] <- wm[t, ] / sum(wm[t, ])
  }
  return(list(thetas = thetas, xs = xs, As = As, wm = wm, Wm = Wm, ess = ess))
}
