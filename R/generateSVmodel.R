
update_Xt <- function(mu, rho, Ut, Xs) {
  Xt <- rho * (Xs - mu) + Ut + mu
  return(Xt)
}

test_generate_SV_data <- function(mu, rho, sigma2, tt) {
  U <- rnorm(N, 0, sqrt(sigma2))
  X0 <- exp(rnorm(1, 0, 1)) # initialise with log-normal rv
  Xt <- rho^tt * X0 + U * rho^c((tt-1):0) + mu * (1-rho^tt)
  return(c(X0,Xt))
}
#' @export generate_SV_data
generate_SV_data  <- function(mu, rho, sigma2, tt) {
  U <- rnorm(N, 0, sqrt(sigma2))
  X0 <- exp(rnorm(1, 0, 1)) # initialise with log-normal rv
  Xt <- c(X0)
  for (t in 2:(tt+1)) {
    Xt <- c(Xt, update_Xt(mu, rho, U[t-1], Xt[t-1]))
  }
  return(Xt)
}

