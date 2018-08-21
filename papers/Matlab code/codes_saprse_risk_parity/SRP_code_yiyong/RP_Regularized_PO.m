function [Result] = RP_Regularized_PO(M, w0, tau, TOParam, CovMatrix, ApproxFun, MaxIter,R,r_c)
% M: a predefined sparse matrix with the i-th row being the same as that of the covariance matrix and 0 elsewhere.
% w0: initial point
% tau: the parameter for the quadratic proximal term 
% TOParam: TOParam.lam1 = lambda_1, TOParam.lam2 = lambda_2
% CovMatrix: covaraince matrix
% ApproxFun: the approximation of the indicator function (log, exp, l_p norm)
% MaxIter: maximal number of iterations
% R: R = sqrtm(CovMatrix)
% r_c: r_c = nu/2*(R\(r_m)) with r_m being the mean return vector

% get trade-off parameters
lam1 = TOParam.lam1;
lam2 = TOParam.lam2;
% nu = TOParam.nu;

% define some auxiliary parameters
NumStocks = length(w0);
eyeNS = eye(NumStocks);
oneNS = ones(NumStocks,1);

% get approximation information
Approx_p = ApproxFun.Approx_p;
Approx_eps = ApproxFun.Approx_eps;
ApproxMethod = ApproxFun.method;

% initialization
gam = 0.9;
w_k = w0;
tmprho = General_Approx(w_k, Approx_p, Approx_eps, ApproxMethod);
tmprho2 = tmprho.^2;
x_k = tmprho2 ./ sum(tmprho2);
tmpSw = CovMatrix*w_k;
tmpwdotSw = w_k.*tmpSw;
theta_k = x_k' * tmpwdotSw;

% set container
rec_val = [];
rec_val = [rec_val PO_Objective(w_k, theta_k, CovMatrix, TOParam, ApproxFun, R, r_c)];

% start the SCA algorithm
count = 0;
while 1
    % We take the form of g_i(w) = w_i(CovMatrix*w)_i
    %% update theta
    % compute theta hat
    tmprho = General_Approx(w_k, Approx_p, Approx_eps, ApproxMethod);
    tmprho2 = tmprho.^2;
    x_k = tmprho2 ./ sum(tmprho2);
    tmpSw = CovMatrix*w_k;
    tmpwdotSw = w_k.*tmpSw;
    theta_hat = x_k' * tmpwdotSw;
    
    %% update w
    tmpVar = w_k'*tmpSw;
    
    tmpM_k = [M{:}]'*w_k;
    tmprho_deri = General_Approx_Deri(w_k, Approx_p, Approx_eps, ApproxMethod);
    
    A_k = reshape(tmpM_k, [NumStocks, NumStocks]) * diag(tmprho) + diag((tmpwdotSw - theta_k) .* tmprho_deri);
    A_k = A_k';
    g_k = (tmpwdotSw - theta_k) .* tmprho;
    
    % get weigths d_k
    d_k = General_Approx_Weights(w_k, Approx_p, Approx_eps, ApproxMethod);
    D_k = diag(d_k);
    
    % Q_k and q_k
    Q_k = lam1 * D_k + lam2 * (A_k'*A_k) + tau * eyeNS;
    q_k = (2*lam2) * A_k'*g_k - 2 * (lam2 * A_k'*A_k + tau * eyeNS) * w_k;
    
    % for Portfolio Mangament problem
    tilde_Q_k = Q_k + R'*R;
    tilde_q_k = q_k - 2*R'*r_c;
    
    % solve sub regularized problem - Mosek
    a = [oneNS'];
    blc = [1];
    buc = [1];
    % blx = [oneNS.*0];
%     blx = [];
    blx = zeros(length(w_k),1);
    bux = [];
    % bux = [oneNS] .* 0.1;
    
    % IMPORTANT: mosek solves QP in the standard form: (1/2)w'*Q*w + w'*q
    [proj] = mskqpopt( (2.*tilde_Q_k),tilde_q_k,a,blc,buc,blx,bux,[],'minimize echo(0)');
    solx = proj.sol.itr.xx;
    
    w_hat = solx(1:NumStocks);
    
    % compute updated theta and w
    % backtracking
    obj_cur = rec_val(end);
    tmpgam = gam;
    while 1
        theta_next = theta_k + tmpgam*(theta_hat - theta_k);
        w_next = w_k + tmpgam*(w_hat - w_k);
        obj_next = PO_Objective(w_next, theta_next, CovMatrix, TOParam, ApproxFun, R, r_c);
        if (obj_next <= obj_cur)
            break;
        else
            tmpgam = tmpgam * 0.9;
        end
    end
    
    rec_val = [rec_val obj_next];
    
    count = count + 1;
    if count > 1
        % diff_val = abs(rec_val(count) -rec_val(count-1))
        if ((norm(w_next - w_k, 2) / norm(w_next, 2) < 10^(-10) && abs(theta_next -theta_k) ./ abs(theta_next) < 1e-10) ...
                || count >= MaxIter)
            break;
        end
    end
    
    % update gam
    gam = gam * (1 - 1e-7*gam);
    w_k = w_next;
    theta_k = theta_next;
end

Result.w = w_k;
Result.rec_val = rec_val;
end