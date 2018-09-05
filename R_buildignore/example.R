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

#res <- riskParityPortfolioCVX(mu, Sigma, tau = 1e-6,
#                              zeta = .1, gamma = .99, ftol = 1e-16, wtol = 1e-10)

res <- riskParityPortfolioGenSolver(Sigma, w0 = rep(1/N, N))
cat(res$init_portfolio_weights, '\n')
cat(res$portfolio_weights, '\n')
cat(sum(res$portfolio_weights))
barplot(t(res$portfolio_weights))
barplot(t(res$risk_contributions))
