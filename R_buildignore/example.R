library(riskParityPortfolio)

Sigma <- diag(c(1:50) ^ 2)
mu <- runif(50)
res <- riskParityPortfolioCVX(mu, Sigma)
print(res)
N_iter <- length(res$obj_fun)
plot(c(1:N_iter), res$obj_fun, type = "b", pch=19, cex=.6,
     col = scales::alpha("black", .5), xlab = "Iteration number",
     ylab = "Objective function")
