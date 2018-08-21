clear all; clc;
load('CovMatrix.mat');  % load the covaraince matrix
% v = randn(20,20);
% CovMatrix = v * v';
NumStocks=size(CovMatrix,2);

% parameter settings
ApproxFun.Approx_p = 0.002;  % For the ApproxFun, please refer to the table of the paper. Note that here the method is fixed to be "Log"
ApproxFun.Approx_eps = 1e-8;
ApproxFun.method = 'Log';  
tau = 100*trace(CovMatrix)/(2*NumStocks);  % the paramter for the proximal quadratic term
w0 = ones(NumStocks,1)/NumStocks;  % initial point
MaxIter = 200;  % maximal number of iterations 

% risk parity portfolio
[w_PRP]=RiskParityPort_Quadratic(w0,CovMatrix,NumStocks,ApproxFun,tau,MaxIter);