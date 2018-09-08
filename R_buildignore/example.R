library(riskParityPortfolio)
library(xts)
library(quantmod)
library(PerformanceAnalytics)

N <- 30
Sigma <- diag((1:N) ^ 2)
Sigma[1, 2] <- 10

w_diag <- 1/sqrt(diag(Sigma))
w_diag <- w_diag/sum(w_diag)

res <- riskParityPortfolioGenSolver(Sigma)
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

