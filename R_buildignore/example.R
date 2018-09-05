library(riskParityPortfolio)
library(xts)
library(quantmod)
library(PerformanceAnalytics)

## set begin-end date and stock namelist
#begin_date <- "2013-01-01"
#end_date <- "2017-08-31"
#stock_namelist <- c("AAPL", "AMD", "ADI",  "ABBV", "AET", "A",  "APD", "AA","CF")
#
## download data from YahooFinance
#prices <- xts()
#for (stock_index in 1:length(stock_namelist))
#  prices <- cbind(prices, Ad(getSymbols(stock_namelist[stock_index], 
#                                        from = begin_date, to = end_date, auto.assign = FALSE)))
## compute log-returns and linear returns
#X_log <- diff(log(prices))[-1]
#X_lin <- (prices/lag(prices) - 1)[-1]
#
## or alternatively...
#X_log <- CalculateReturns(prices, "log")[-1]
#X_lin <- CalculateReturns(prices)[-1]
#
#N <- ncol(X_log)  # number of stocks
#T <- nrow(X_log)  # number of days
#
## split data into training and set data
#T_trn <- round(0.7*T)  # 70% of data
#X_log_trn <- X_log[1:T_trn, ]
#X_log_tst <- X_log[(T_trn+1):T, ]
#X_lin_trn <- X_lin[1:T_trn, ]
#X_lin_tst <- X_lin[(T_trn+1):T, ]
#
#mu <- colMeans(X_log_trn)
#Sigma <- cov(X_log_trn)

N <- 10
mu <- runif(N)
Sigma <- diag((1:N) ^ 2)
Sigma[1, 2] <- 50
#res <- riskParityPortfolioCVX(mu, Sigma, tau = 1e-6,
#                              zeta = .1, gamma = .99, ftol = 1e-16, wtol = 1e-10)

w_diag <- 1/sqrt(diag(Sigma))
w_diag <- w_diag/sum(w_diag)
start <- Sys.time()
res <- riskParityPortfolioGenSolver(Sigma)
end <- Sys.time()
print(end - start)
start <- Sys.time()
res_cvx <- riskParityPortfolioCVX(Sigma)
end <- Sys.time()
print(end - start)
start <- Sys.time()
res_qp <- riskParityPortfolioQP(Sigma)
end <- Sys.time()
print(end - start)
cat(res$init_portfolio_weights, '\n')
cat(res$portfolio_weights, '\n')
cat(res_cvx$portfolio_weights, '\n')
cat(w_diag, '\n')
cat(sum(res$portfolio_weights))
barplot(t(res$portfolio_weights))
barplot(t(res$risk_contributions))
barplot(t(res_cvx$portfolio_weights))
barplot(t(res_cvx$risk_contributions))
barplot(t(res_qp$portfolio_weights))
barplot(t(res_qp$risk_contributions))
barplot(t(w_diag))
barplot(t(w_diag * (Sigma %*% w_diag)))
Niter <- length(res_qp$obj_fun)
plot(c(1:Niter), res_qp$obj_fun)

