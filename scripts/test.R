detach("package:smc", unload=TRUE)
rm(list = ls())
library(smc)
library(Rcpp)

set.seed(1)

## Bootstrap filter
tmax <- 1000
mu <- -1
rho <- 0.95
sigma <- 0.15

Xt <- generate_SV_data(mu, rho, sigma, tmax)
Yt <- as.matrix(rnorm(tmax+1, mean = 0, sd = sqrt(exp(Xt))))

boot_sv <- Bootstrap_SV$new(data = Yt, mu = -1, sigma = 0.15, rho = 0.95)

N <- 5000
output <- bootstrap_filter(boot_sv, N, tmax)

mod <- Module("particles", PACKAGE="smc")
Bootstrap_SV_C <- mod$Bootstrap_SV_C
boot_sv_rcpp <- new(Bootstrap_SV_C, Yt, -1, 0.15, 0.95)
output2 <- mod$bootstrap_filter_rcpp(boot_sv_rcpp, N ,tmax)

par(mfrow = c(1,2))
plot(Xt, type = "l")
lines(output$mx, col = "red")
lines(output2$mx, col="green")
plot(1:tmax, output$ess, type = "l")

rm(output)
rm(output2)
gc()

library(microbenchmark)

run_filter <- function(filter_fn, ...) {
  out <- filter_fn(...)
  rm(out)
  gc()
  return(NA)
}
funcs <- c(call("run_filter", bootstrap_filter, boot_sv, N, tmax), call("run_filter", mod$bootstrap_filter_rcpp, boot_sv_rcpp, N, tmax))

bench_res <- microbenchmark(list=setNames(funcs, c("R","Rcpp")),
               times=5L)

print(bench_res, signif=3)


## SMC^2
tmax <- 1000
mu <- -1
rho <- 0.95
sigma <- 0.15
Xt <- generate_SV_data(mu, rho, sigma, tmax)
Yt <- as.matrix(rnorm(tmax+1, mean = 0, sd = sqrt(exp(Xt))))
Nx <- 100
Nt <- 1000
sd_prior <- 0.2
mu_prior <- -1
mu_prior <- -0.7
sd_prop <- 1.5

smc_results <- smc_squared(Yt, Nx, Nt, sigma, rho, 
                           mu_prior = -0.7, sd_prior = 0.2, sd_prop = 1)
                           mu_prior = mu_prior, sd_prior = sd_prior, sd_prop = sd_prop)

## Get mean and sd of theta from last iteration and then run a bootstrap filter

mut = sum(smc_results$Wm[tmax+1,] * smc_results$thetas[tmax+1, ])/sum(smc_results$Wm[tmax+1,])
sd_t = sqrt(sum(smc_results$Wm[tmax+1,] * (smc_results$thetas[tmax+1, ] - mut)^2)/sum(smc_results$Wm[tmax+1,]))

# smc_results <- smc_squared(Yt, Nx, Nt, sigma, rho, 
#                            mu_prior = mut, sd_prior = 0.2, sd_prop = sd_t)
boot_sv_rcpp <- new(Bootstrap_SV_C, Yt, mut, 0.15, 0.95)
output2 <- mod$bootstrap_filter_rcpp(boot_sv_rcpp, N ,tmax)
