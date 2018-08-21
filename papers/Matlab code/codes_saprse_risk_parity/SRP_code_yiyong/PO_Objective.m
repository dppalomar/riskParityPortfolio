function [Fval] = PO_Objective(w_k, theta_k, CovMatrix, TOParam, ApproxFun, R, r_c)

% trade-off parameters
lam1 = TOParam.lam1;
lam2 = TOParam.lam2;

Approx_p = ApproxFun.Approx_p;
Approx_eps = ApproxFun.Approx_eps;
ApproxMethod = ApproxFun.method;

tmpSw = CovMatrix*w_k;
tmpwdotSw = w_k.*tmpSw;


tmprho = General_Approx(w_k, Approx_p, Approx_eps, ApproxMethod);
% Fterm = - nu*w_k'*MeanVec + tmpVar;
Fterm = norm(r_c-R*w_k)^2;
Secterm = sum(tmprho);
tmpSterm = (tmpwdotSw - theta_k) .* tmprho;
Sterm = sum(tmpSterm.^2);

% Fval = Fterm + lam * Sterm;
Fval = Fterm + lam1 * Secterm + lam2 * Sterm;
end