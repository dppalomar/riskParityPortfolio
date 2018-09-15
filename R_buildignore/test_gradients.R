library(numDeriv)  #install.packages("numDeriv")
#set.seed(123)


# generate a random Sigma and w
N <- 5
Sigma12 <- matrix(rnorm(N^2), N, N)
Sigma <- Sigma12 %*% t(Sigma12)
w <- runif(N)

# basic quantities
b <- rep(1/N, N)
Sigma_w <- as.vector(Sigma %*% w)
r <- w*Sigma_w
sum_r <- sum(r)



# Formulation “rc double-index”
risk <- function(w) {
  Sigma_w <- as.vector(Sigma %*% w)
  r <- w*Sigma_w
  2*N*sum(r^2) - 2*sum(r)^2
}
v <- N*r-sum(r)
risk_grad <- as.vector(4*(Sigma %*% (w*v) + Sigma_w*v))
risk_grad_num <- grad(risk, w)
norm(risk_grad - risk_grad, "2")

g <- function(w) {
  N <- length(w)
  Sigma_w <- as.vector(Sigma %*% w)
  r <- w*Sigma_w
  rep(r, times = N) - rep(r, each = N)
}
Ut <- diag(Sigma_w) + Sigma*w
g_jac <- matrix(rep(t(Ut), N), ncol = N, byrow = TRUE) - matrix(rep(Ut, each = N), ncol = N)
g_jac_num <- jacobian(g, w)
norm(g_jac - g_jac_num, "F")



# Formulation “rc-over-var vs b”
risk <- function(w) {
  Sigma_w <- as.vector(Sigma %*% w)
  r <- w*Sigma_w
  sum((r/sum(r) - b)^2)
}
r_b <- r/sum_r - b
v <- r_b - as.numeric(r_b %*% r)/sum_r
risk_grad <- (2/sum_r) * as.vector(Sigma %*% (w*v) + Sigma_w*v)
risk_grad_num <- grad(risk, w)
norm(risk_grad - risk_grad, "2")

g <- function(w) {
  Sigma_w <- as.vector(Sigma %*% w)
  r <- w*Sigma_w
  r/sum(r) - b
}
Ut <- diag(Sigma_w) + Sigma*w
g_jac <- Ut/sum_r - 2/(sum_r^2) * r %o% Sigma_w
g_jac_num <- jacobian(g, w)
norm(g_jac - g_jac_num, "F")

v <- colSums(Ut)
g_jac_bis <- Ut/sum_r - 1/(sum_r^2) * r %o% v
norm(g_jac_bis - g_jac_num, "F")



# Formulation "rc-over-sd vs b-times-sd"
risk <- function(w) {
  Sigma_w <- as.vector(Sigma %*% w)
  r <- w*Sigma_w
  sqrt_sum_r <- sqrt(sum(r))
  sum((r/sqrt_sum_r - b*sqrt_sum_r)^2)
}
r_b <- r/sum_r - b
v <- r_b - as.numeric(r_b %*% r)/sum_r
risk_grad <- as.vector(Sigma %*% (w*v) + Sigma_w*v)
risk_grad_num <- grad(risk, w)
norm(risk_grad - risk_grad, "2")

g <- function(w) {
  Sigma_w <- as.vector(Sigma %*% w)
  r <- w*Sigma_w
  sqrt_sum_r <- sqrt(sum(r))
  r/sqrt_sum_r - b*sqrt_sum_r
}
Ut <- diag(Sigma_w) + Sigma*w
g_jac <- Ut/sqrt(sum_r) - (r/(sum_r^(3/2)) + b/sqrt(sum_r)) %o% Sigma_w
g_jac_num <- jacobian(g, w)
norm(g_jac - g_jac_num, "F")



g <- function(r) r/sum(r)
jacobian(g, r)
(1/sum_r) * (diag(N) - (1/sum_r)*matrix(t(r), N, N, byrow = TRUE))

g1 <- function(r) r[1]/sum(r)
grad(g1, r)
(1/sum_r) * (c(1, 0, 0, 0, 0) - r/sum_r)



x <- c(1, 2)
f <- function(x) x/sum(x)
jacobian(f, x)
(1/sum(x)) * (diag(2) - (1/sum(x))*matrix(x, 2, 2))
