detach("package:smc", unload=TRUE)
rm(list = ls())
library(smc)
library(Rcpp)

set.seed(1)

tmax <- 500
mu <- -1
rho <- 0.95
sigma <- 0.15

Xt <- generate_SV_data(mu, rho, sigma, tmax)
Yt <- as.matrix(rnorm(tmax, mean = 0, sd = sqrt(exp(Xt))))
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

bench_res <- microbenchmark(rv <- run_filter(bootstrap_filter, boot_sv, N, tmax),
               cv <- run_filter(mod$bootstrap_filter_rcpp, boot_sv_rcpp, N, tmax),
               times=5L)

