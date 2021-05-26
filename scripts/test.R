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
set.seed(2)
output <- bootstrap_filter(boot_sv, N, tmax)

boot_sv_rcpp <- new(Bootstrap_SV_C, Yt, -1, 0.15, 0.95)
set.seed(2)
output2 <- bootstrap_filter_rcpp(boot_sv_rcpp, N ,tmax)

par(mfrow = c(2, 1))
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
Nt <- 200
sd_prior <- 0.2
mu_prior <- -1.1
sd_prop = 0.5

smc_results <- smc_squared_rcpp(Yt, Nx, Nt, sigma, rho, mu_prior, 
                                sd_prior, sd_prop)
