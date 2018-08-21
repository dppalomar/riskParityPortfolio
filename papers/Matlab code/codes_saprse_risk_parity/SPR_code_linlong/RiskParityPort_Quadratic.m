function [w_PRP, Lam1sQ_best, Alpha2sQ_best]=RiskParityPort_Quadratic(Data_intest,CovMatrix_train,MeanVec_train,CovMatrix,MeanVec,Lam1s,Alpha2s,nu,NumStocks)
% Data_intest: the daily return matrix 
% CovMatrix_train: the covariance matrix for training data
% MeanVec_train: the mean vector for training data
% CovMatrix: the covariance matrix for the covariance matrix for Data_intest
% MeanVec: the  mean vector for Data_intest
% Lam1s: lambda_1 in the problem formulation
% Alpha2s: lambda_2 in the problem formualtion
% nu: the parameter for the mean-variance goal, w'*\Sigma*w-nu*w'*mu
% NumStocks: the number of stockes
%% approximation settings
oneN = ones(NumStocks, 1);
ApproxFun.Approx_p = 0.002;
ApproxFun.Approx_eps = 1e-8;
ApproxFun.method = 'Log';
TOParam.lam1 = 0;
TOParam.lam2 = 0;
TOParam.nu = 0;

%% cross validation to choose the paraters lambda1 and lambda2
R=sqrtm(CovMatrix_train);
r_index=nu/2*(R\(MeanVec_train'));
NumLam1 = length(Lam1s);
NumAlpha2 = length(Alpha2s);
Lam1sQ_best=[];
Alpha2sQ_best=[];
SR_PRP=[];
PortTemp=[];
for whichlam1 = 1:NumLam1
    TOParam.lam1 = Lam1s(whichlam1);
    for whichlam2 = 1:NumAlpha2
        TOParam.lam2 = Alpha2s(whichlam2);
        ScaleCovMatrix = 1e4 .* CovMatrix_train;
        M = cell(NumStocks,1);
        for whichi = 1:NumStocks
            tmpM = sparse(NumStocks,NumStocks);
            tmpM(whichi,:) = ScaleCovMatrix(whichi,:);
            M{whichi} = (tmpM + tmpM');
        end
        tau = 0.01*trace(ScaleCovMatrix)/(2*NumStocks);
        w0 = oneN ./ NumStocks;
        tmp = RP_Regularized_Quadratic(M, w0, tau, TOParam, ScaleCovMatrix, ApproxFun, 200,R,r_index);
        % calculate the SR and update Lam1s_bestQ & Alpha2s_bestQ
        PortTemp.PRP{1}=tmp.w;   
        for k=1:size(Data_intest,1)
            temp=(1+Data_intest(k,:))'.*PortTemp.PRP{k};
            PortTemp.PRP{k+1}=temp./sum(temp);
        end  
        re_PRP=[];
        for k=1:size(Data_intest,1)
            re_PRP=[re_PRP Data_intest(k,:)*PortTemp.PRP{k}];
        end
        return_PRP = sum(re_PRP)/length(re_PRP);
        vol_PRP = norm(re_PRP-return_PRP)/sqrt(length(re_PRP)-1);
        %             SR_PRP=[SR_PRP return_PRP/vol_PRP];
        %             SR_PRP=[SR_PRP nu*return_PRP-(vol_PRP)^2];
        money=cumprod(1+re_PRP);
        SR_PRP=[SR_PRP 0.05*money(end)+(nu*return_PRP-(vol_PRP)^2)];
        if isempty(Lam1sQ_best)
            Lam1sQ_best=TOParam.lam1;
            Alpha2sQ_best=TOParam.lam2;
        else if SR_PRP(end)>SR_PRP(end-1)
                Lam1sQ_best=TOParam.lam1;
                Alpha2sQ_best=TOParam.lam2;
            end
        end
    end
end

%% design the portfolio weights
R=sqrtm(CovMatrix);
r_index=nu/2*(R\(MeanVec'));
TOParam.lam1=Lam1sQ_best;
TOParam.lam2=Alpha2sQ_best;
ScaleCovMatrix = 1e4 .* CovMatrix;
M = cell(NumStocks,1);
for whichi = 1:NumStocks
    tmpM = sparse(NumStocks,NumStocks);
    tmpM(whichi,:) = ScaleCovMatrix(whichi,:);
    M{whichi} = (tmpM + tmpM');
end
tau = 0.01*trace(ScaleCovMatrix)/(2*NumStocks);
w0 = oneN ./ NumStocks;
tmp = RP_Regularized_Quadratic(M, w0, tau, TOParam, ScaleCovMatrix, ApproxFun, 200,R,r_index);
w_PRP=tmp.w;