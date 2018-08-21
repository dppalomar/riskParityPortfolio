function [w_PRP]=RiskParityPort_Quadratic(w0,CovMatrix,NumStocks,ApproxFun,tau, MaxIter)
TOParam.lam1=0;
TOParam.lam2=1;

ScaleCovMatrix = 1e4 .* CovMatrix;
M = cell(NumStocks,1);
for whichi = 1:NumStocks
    tmpM = sparse(NumStocks,NumStocks);
    tmpM(whichi,:) = ScaleCovMatrix(whichi,:);
    M{whichi} = (tmpM + tmpM');
end
% tau = 0.01*trace(ScaleCovMatrix)/(2*NumStocks);
% w0 = oneN ./ NumStocks;
w_PRP = RP_Regularized_Quadratic(M, w0, tau, TOParam, ScaleCovMatrix, ApproxFun, MaxIter);
end