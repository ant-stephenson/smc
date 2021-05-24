detach("package:smc", unload=TRUE)
rm(list = ls())
library(smc)

set.seed(1)

# # not the same for whatever reason, so use test version
# test_Xt <- test_generate_SV_data(mu, rho, sigma2, tt)
# Xt <- generate_SV_data(mu, rho, sigma2, tt)
# Yt <- rnorm(tt, 0, 1) * exp(Xt)
# 
# print(test_Xt)
# print(Xt)

tmax <- 100
mu <- -1
rho <- 0.95
sigma <- 0.15

Xt <- generate_SV_data(mu, rho, sigma, tmax)
Yt <- as.matrix(rnorm(tmax+1, mean = 0, sd = sqrt(exp(Xt))))

boot_sv <- Bootstrap_SV$new(data = Yt, mu = -1, sigma = 0.15, rho = 0.95)

N <- 10
output <- bootstrap_filter(boot_sv, N, tmax)

par(mfrow = c(2, 1))
#plot(Yt, type = "l")
plot(Xt, type = "l")
lines(output$mx, col = "red")
#plot(output$sdx, type = "l")
plot(1:tmax, output$ess, type = "l")
rm(output)

tmax <- 1000
mu <- -1
rho <- 0.95
sigma <- 0.15
Xt <- generate_SV_data(mu, rho, sigma, tmax)
Yt <- as.matrix(rnorm(tmax+1, mean = 0, sd = sqrt(exp(Xt))))
Nx <- 100
Nt <- 1000
sd_prior <- 0.2
mu_prior <- 0

smc_results <- smc_squared(Yt, Nx, Nt, sigma, rho, 
                           mu_prior = -1, sd_prior = 1, sd_prop = 1)
