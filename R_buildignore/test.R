library(xts)
library(riskParityPortfolio)

# data
Sigma <- rbind(c(1.0000, 0.0015, -0.0119),
               c(0.0015, 1.0000, -0.0308),
               c(-0.0119, -0.0308, 1.0000))
b <- c(0.1594, 0.0126, 0.8280)
N <- 3


# Formulation "rc-double-index"
res_gen_solver1 <- riskParityPortfolioGenSolver(Sigma, formulation = "rc-double-index")
res_gen_solver1$w
last(res_gen_solver1$obj_fun)

res_sca1 <- riskParityPortfolioSCA(Sigma, formulation = "rc-double-index")
res_sca1$w
last(res_sca1$obj_fun)

# Formulation "rc-over-var-vs-b"
res_gen_solver2 <- riskParityPortfolioGenSolver(Sigma, b, formulation = "rc-over-var-vs-b")
res_gen_solver2$w
last(res_gen_solver2$obj_fun)

res_sca2 <- riskParityPortfolioSCA(Sigma, b, formulation = "rc-over-var-vs-b")
res_sca2$w
last(res_sca2$obj_fun)

# Formulation "rc-over-sd-vs-b-times-sd"
res_gen_solver3 <- riskParityPortfolioGenSolver(Sigma, b, formulation = "rc-over-sd-vs-b-times-sd")
res_gen_solver3$w
last(res_gen_solver3$obj_fun)

res_sca3 <- riskParityPortfolioSCA(Sigma, b, formulation = "rc-over-sd-vs-b-times-sd")
res_sca3$w
last(res_sca3$obj_fun)

