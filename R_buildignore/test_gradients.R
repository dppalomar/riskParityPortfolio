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


#
#  Formulation “rc-double-index”
#
R_rc_double_index(w, Sigma)
risk_grad <- R_grad_rc_double_index(w, Sigma)
risk_grad_num <- grad(R_rc_double_index, x = w, Sigma = Sigma)
norm(risk_grad - risk_grad, "2")

g_jac <- A_rc_double_index(w, Sigma)
g_jac_num <- jacobian(g_rc_double_index, x = w, Sigma = Sigma)
norm(g_jac - g_jac_num, "F")


#
#  Formulation “rc-over-b-double-index”
#
R_rc_over_b_double_index(w, Sigma, b)
risk_grad <- R_grad_rc_over_b_double_index(w, Sigma, b)
risk_grad_num <- grad(R_rc_over_b_double_index, x = w, Sigma = Sigma, b = b)
norm(risk_grad - risk_grad, "2")

g_jac <- A_rc_over_b_double_index(w, Sigma, b)
g_jac_num <- jacobian(g_rc_over_b_double_index, x = w, Sigma = Sigma, b = b)
norm(g_jac - g_jac_num, "F")


#
#  Formulation “rc-over-var-vs-b”
#
R_rc_over_var_vs_b(w, Sigma, b)
risk_grad <- R_grad_rc_over_var_vs_b(w, Sigma, b)
risk_grad_num <- grad(R_rc_over_var_vs_b, x = w, Sigma = Sigma, b = b)
norm(risk_grad - risk_grad, "2")

g_jac <- A_rc_over_var_vs_b(w, Sigma, b)
g_jac_num <- jacobian(g_rc_over_var_vs_b, x = w, Sigma = Sigma, b = b)
norm(g_jac - g_jac_num, "F")


#
#  Formulation “rc-over-var”
#
R_rc_over_var(w, Sigma)
risk_grad <- R_grad_rc_over_var(w, Sigma)
risk_grad_num <- grad(R_rc_over_var, x = w, Sigma = Sigma)
norm(risk_grad - risk_grad, "2")

g_jac <- A_rc_over_var(w, Sigma)
g_jac_num <- jacobian(g_rc_over_var, x = w, Sigma = Sigma)
norm(g_jac - g_jac_num, "F")


#
#  Formulation "rc-over-sd vs b-times-sd"
#
R_rc_over_sd_vs_b_times_sd(w, Sigma, b)
risk_grad <- R_grad_rc_over_sd_vs_b_times_sd(w, Sigma, b)
risk_grad_num <- grad(R_rc_over_sd_vs_b_times_sd, x = w, Sigma = Sigma, b = b)
norm(risk_grad - risk_grad, "2")

g_jac <- A_rc_over_sd_vs_b_times_sd(w, Sigma, b)
g_jac_num <- jacobian(g_rc_over_sd_vs_b_times_sd, x = w, Sigma = Sigma, b = b)
norm(g_jac - g_jac_num, "F")


#
#  Formulation "rc vs b-times-var"
#
R_rc_vs_b_times_var(w, Sigma, b)
risk_grad <- R_grad_rc_vs_b_times_var(w, Sigma, b)
risk_grad_num <- grad(R_rc_vs_b_times_var, x = w, Sigma = Sigma, b = b)
norm(risk_grad - risk_grad, "2")

g_jac <- A_rc_vs_b_times_var(w, Sigma, b)
g_jac_num <- jacobian(g_rc_vs_b_times_var, x = w, Sigma = Sigma, b = b)
norm(g_jac - g_jac_num, "F")


#
#  Formulation "rc vs theta"
#
theta <- sum(r)/N  #BTW, I would use this as the initial value
R_rc_vs_theta(c(w, theta), Sigma)
risk_grad <- R_rc_vs_theta(c(w, theta), Sigma)
risk_grad_num <- grad(R_rc_vs_theta, x = c(w, theta), Sigma = Sigma)
norm(risk_grad - risk_grad, "2")

g_jac <- A_rc_vs_theta(c(w, theta), Sigma)
g_jac_num <- jacobian(g_rc_vs_theta, x = c(w, theta), Sigma = Sigma)
norm(g_jac - g_jac_num, "F")

