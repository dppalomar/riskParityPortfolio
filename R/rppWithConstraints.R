# this implementation is not optimized, it's just a first attempt
rpp_equality_constraints_iteration_R <- function(Cmat, cvec, Qk, qk) {
  inv_Qk <- solve(Qk)  # it's faster to compute the inverse!
  lmd_k <- -solve(Cmat %*% inv_Qk %*% t(Cmat), Cmat %*% inv_Qk %*%qk + cvec)
  #lmd_k <- -solve(Cmat %*% solve(Qk, t(Cmat)), Cmat %*% solve(Qk, qk) + cvec)
  w_hat <- -inv_Qk %*% (qk + t(Cmat) %*% lmd_k)
  #w_hat <- -solve(Qk, qk + t(Cmat) %*% lmd_k)
  return(as.vector(w_hat))
}


rpp_eq_and_ineq_constraints_iteration_R <- function(Cmat, cvec, Dmat, dvec, Qk,
                                                    qk, wk,
                                                    dual_mu_0, dual_mu_minus_1,
                                                    dual_lmd_0, dual_lmd_minus_1,
                                                    wtol = 1e-6, maxiter = 500) {
  B <- rbind(Cmat, Dmat)
  inv_Qk <- solve(Qk)
  L <- norm(B %*% inv_Qk %*% t(B), type = "2")
  #L <- norm(B %*% solve(Qk, t(B)), type = "2")
  w_tilde_i_minus_1 <- wk
  dual_lmd_i <- dual_lmd_0;                dual_mu_i <- dual_mu_0
  dual_lmd_i_minus_1 <- dual_lmd_minus_1;  dual_mu_i_minus_1 <- dual_mu_minus_1
  for (i in c(0:maxiter)) {
    w_tilde_i <- -inv_Qk %*% (qk + t(Cmat) %*% dual_lmd_i + t(Dmat) %*% dual_mu_i)
    #w_tilde_i <- -solve(Qk, qk + t(Cmat) %*% dual_lmd_i + t(Dmat) %*% dual_mu_i)
    w_tilde_bar_i <- w_tilde_i + (i-1)/(i+2) * (w_tilde_i - w_tilde_i_minus_1)
    dual_lmd_i_plus_1 <- dual_lmd_i + (i-1)/(i+2) * (dual_lmd_i - dual_lmd_i_minus_1) + 1/L*(Cmat %*% w_tilde_bar_i - cvec)
    dual_mu_i_plus_1  <- pmax(0, dual_mu_i  + (i-1)/(i+2) * (dual_mu_i  - dual_mu_i_minus_1 ) + 1/L*(Dmat %*% w_tilde_bar_i - dvec))
    if (i > 0 && all(abs(w_tilde_bar_i - w_tilde_bar_i_minus_1) <=
                     .5 * wtol * (abs(w_tilde_bar_i) + abs(w_tilde_bar_i_minus_1)))) break;
    w_tilde_i_minus_1 <- w_tilde_i
    w_tilde_bar_i_minus_1 <- w_tilde_bar_i
    dual_lmd_i_minus_1 <- dual_lmd_i;  dual_mu_i_minus_1 <- dual_mu_i
    dual_lmd_i <- dual_lmd_i_plus_1;   dual_mu_i <- dual_mu_i_plus_1;
  }
  return(list(dual_mu_i_plus_1, dual_mu_i, dual_lmd_i_plus_1, dual_lmd_i, w_tilde_bar_i))
}
