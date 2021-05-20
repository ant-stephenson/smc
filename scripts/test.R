detach("package:smc", unload=TRUE)
rm(list = ls())
library(smc)

set.seed(1)

tmax <- 500
mu <- -1
rho <- 0.95
sigma <- 0.15

Xt <- generate_SV_data(mu, rho, sigma, tmax)
Yt <- as.matrix(rnorm(tmax, mean = 0, sd = sqrt(exp(Xt))))
boot_sv <- Bootstrap_SV$new(data = Yt, mu = -1, sigma = 0.15, rho = 0.95)

N <- 50000
output <- bootstrap_filter(boot_sv, N, tmax)

par(mfrow = c(1,2))
plot(Xt, type = "l")
lines(output$mx, col = "red")
plot(1:500, output$ess, type = "l")
rm(output)
