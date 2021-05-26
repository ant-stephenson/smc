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

# library(microbenchmark)
# 
# run_filter <- function(filter_fn, ...) {
#   out <- filter_fn(...)
#   rm(out)
#   gc()
#   return(NA)
# }
# funcs <- c(call("run_filter", bootstrap_filter, boot_sv, N, tmax), call("run_filter", mod$bootstrap_filter_rcpp, boot_sv_rcpp, N, tmax))
# 
# bench_res <- microbenchmark(list=setNames(funcs, c("R","Rcpp")),
#                times=5L)
# 
# print(bench_res, signif=3)


## SMC^2 -run once for burn-in/parameters
tmax <- 1000
mu <- -1
rho <- 0.95
sigma <- 0.15
Xt <- generate_SV_data(mu, rho, sigma, tmax)
Yt <- as.matrix(rnorm(tmax+1, mean = 0, sd = sqrt(exp(Xt))))
Nx <- 1000
Nt <- 100
sd_prior <- 0.2
mu_prior <- -0.7
sd_prop <- 1.5

mu_prior = c(-0.7, 0.2)
sd_prior <- c(0.2, 0.05)
sd_prop <- c(1.5, 1.0)

smc_results <- smc_squared(Yt, Nx, Nt, sigma, rho, 
                           mu_prior = mu_prior, sd_prior = sd_prior, sd_prop = sd_prop)

# check convergence of parameters
plot(rowMeans(smc_results$thetas[,,1]), type="l")
lines(rep(mu, tmax+1), col="red")
# miny <- min(sigma, min(rowMeans(smc_results$thetas[,,2])))
# maxy <- max(sigma, max(rowMeans(smc_results$thetas[,,2])))
plot(rowMeans(smc_results$thetas[,,2]), type="l", ylim=c(miny, maxy), log="y")
lines(rep(sigma, tmax+1), col="red")

## Get mean and sd of theta from last iteration and then run a bootstrap filter

mut = sum(smc_results$Wm[tmax+1,] * smc_results$thetas[tmax+1, ])/sum(smc_results$Wm[tmax+1,])
sd_t = sqrt(sum(smc_results$Wm[tmax+1,] * (smc_results$thetas[tmax+1, ] - mut)^2)/sum(smc_results$Wm[tmax+1,]))

# smc_results <- smc_squared(Yt, Nx, Nt, sigma, rho, 
#                            mu_prior = mut, sd_prior = 0.2, sd_prop = sd_t)
boot_sv_rcpp <- new(Bootstrap_SV_C, Yt, mut, 0.15, 0.95)
output2 <- mod$bootstrap_filter_rcpp(boot_sv_rcpp, N ,tmax)
