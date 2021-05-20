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

tmax <- 500
mu <- -1
rho <- 0.95
sigma <- 0.15

Xt <- generate_SV_data(mu, rho, sigma, tmax)
Yt <- as.matrix(rnorm(tmax+1, 0, 1) * exp(Xt))
boot_sv <- Bootstrap_SV$new(data = Yt, mu = -0, sigma = 1, rho = 0.8)

N <- 1000
output <- bootstrap_filter(boot_sv, N, tmax)
str(output)

plot(Yt, type = "l")
lines(output$mx, col = "red")
