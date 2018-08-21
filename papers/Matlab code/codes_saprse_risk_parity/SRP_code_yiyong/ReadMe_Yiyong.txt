The two functions are RP_Regularized_PO.m and RP_Regularized_PO_IterativeL1.m.
The parameter settings are as follows (use RP_Regularized_PO.m as the example):

RP_Regularized_PO.m:
Input:
% M: a predefined sparse matrix with the i-th row being the same as that of the covariance matrix and 0 elsewhere.
% w0: initial point
% tau: the parameter for the quadratic proximal term 
% TOParam: TOParam.lam1 = lambda_1, TOParam.lam2 = lambda_2
% CovMatrix: covaraince matrix
% ApproxFun: the approximation of the indicator function (log, exp, l_p norm)
% MaxIter: maximal number of iterations
% R: R = sqrtm(CovMatrix)
% r_c: r_c = nu/2*(R\(r_m)) with r_m being the mean return vector
Output:
% Result: where Result.w is the designed portfolio.

The parameters need to be set are TOParam and tau.

Since the shared folder between Yiyong and I has disappeared, I don't have the sample test or the main fucntion for RP_Regularized_PO.m. Yiyong has the original codes.