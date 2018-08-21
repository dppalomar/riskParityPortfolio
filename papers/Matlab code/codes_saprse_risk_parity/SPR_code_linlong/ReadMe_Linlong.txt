sample_main.m is the sample test file. It bascially cantains all the necessary parts (dived the data into rolliing windows, the benchmarks, obtain the final portfolio vector, .etc). In the sample_main, you only need to focus on line 48 - 60 corresponding to the risk parity algorithms.

Because we derived to approxiamtions, one is quadratic, the other is linear. In the following, I will use the quadratic one as the example.

First, RiskParityPort_Quadratic.m is a relatively compelete part, only requeiring the input data and some parameters, then it will output the portfolio vector. The meaning of the input parameters are explained in line 2-10 of RiskParityPort_Quadratic.m. Lam1s, Alpha2s and nu are the parameters need to set by the user. In addition, we call the function RP_Regualrized_Quadratic.m in line 45 of RiskParityPort_Quadratic.m. In RP_Regualrized_Quadratic.m, one parameter named gam in line 18 need to be set. The realted codes are line 18, 85-95 and 109.The other parameter is tau, which is the parameter for the proximal term in Aldo's framework. Tau is already set as input parameter of RiskParityPort_Quadratic.m and the related line is line59-60.

Second, in RiskParity_Quadratic.m, line 23-68 is the code of cross validation for lambda_1 and lambda_2.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Summary: 

RiskParityPort_Quadratic.m:

%%%  Input  %%%%
Data_intest: the daily return matrix 
CovMatrix_train: the covariance matrix for training data
MeanVec_train: the mean vector for training data
CovMatrix: the covariance matrix for the covariance matrix for Data_intest
MeanVec: the  mean vector for Data_intest
Lam1s: lambda_1 in the problem formulation
Alpha2s: lambda_2 in the problem formualtion
nu: the parameter for the mean-variance goal, w'*\Sigma*w-nu*w'*mu
NumStocks: the number of stockes

%%%  Output  %%%%
portfolio vector for sparse risk parity


Five parameters need to be set/tuned:
Related to the problem formulation od spare risk parity:
Lam1s, Alpha2s, nu
Related to the Aldo's framework:
tau, gam
