update_Xt <- function(mu, rho, Ut, Xs) {
  Xt <- rho * (Xs - mu) + Ut + mu
  return(Xt)
}

test_generate_SV_data <- function(mu, rho, sigma, tmax) {
  U <- rnorm(tmax, 0, sqrt(sigma))
  X0 <- exp(rnorm(1, 0, 1)) 
  Xt <- rho^tmax * X0 + U * rho^c((tmax-1):0) + mu * (1-rho^tmax)
  return(c(X0,Xt))
}
#' @export generate_SV_data
generate_SV_data  <- function(mu, rho, sigma, tmax) {
  U <- rnorm(tmax, 0, sigma)
  X0 <- rnorm(1, 0, 1) 
  Xt <- c(X0)
  for (t in 2:(tmax+1)) {
    Xt <- c(Xt, update_Xt(mu, rho, U[t-1], Xt[t-1]))
  }
  return(Xt)
}

