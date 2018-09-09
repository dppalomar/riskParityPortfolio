library(riskParityPortfolio)
library(xts)
library(quantmod)
library(PerformanceAnalytics)

N <- 50
Sigma <- diag((1:N) ^ 2)
corr <- matrix(runif(N ^2), nrow = N)
Sigma[1, 2] <- 10
Sigma[2, 1] <- 10
Sigma <- Sigma + corr + t(corr)

w_diag <- 1/sqrt(diag(Sigma))
w_diag <- w_diag/sum(w_diag)

res <- riskParityPortfolioGenSolver(Sigma, ftol = 1e-9, wtol = 1e-9)
res_qp <- riskParityPortfolioQP(Sigma)
w_all <- c("risk-parity-gen-solver" = res$w,
           "risk-parity-qp" = res_qp$w,
           "risk-parity-diag" = w_diag)
barplot(t(res$w))
barplot(t(res$risk_contributions))
barplot(t(res_qp$w))
barplot(t(res_qp$risk_contributions))
barplot(t(w_diag))
barplot(t(w_diag * (Sigma %*% w_diag)))

plot(res_qp$elapsed_time, res_qp$obj_fun)
plot(res$elapsed_time, res$obj_fun)

