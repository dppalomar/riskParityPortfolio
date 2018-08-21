function [w_PRP_L1, Lam1sL_best, Alpha2sL_best]=RiskParityPort_Linear(Data_intest,CovMatrix_train,MeanVec_train,CovMatrix,MeanVec,Lam1s,Alpha2s,nu,NumStocks)
%% parameter settings for approxiamtion
oneN = ones(NumStocks, 1);
ApproxFun.Approx_p = 0.002;
ApproxFun.Approx_eps = 1e-15;
ApproxFun.method = 'Log';
TOParam.lam1 = 0;
TOParam.lam2 = 0;
TOParam.nu = 0;

%% cross validation 
R=sqrtm(CovMatrix_train);
r_index=nu/2*(R\(MeanVec_train'));
NumLam1 = length(Lam1s);
NumAlpha2 = length(Alpha2s);
Lam1sL_best=[];
Alpha2sL_best=[];
SR_PRP_L1=[];
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
        tmp = RP_Regularized_Linear(M, w0, tau, TOParam, ScaleCovMatrix, ApproxFun, 200, R, r_index);
        % calculate the SR of Data_test and update Lam1s_bestQ & Alpha2s_bestQ
        PortTemp.PRP_L1{1}=tmp.w; 
        for k=1:size(Data_intest,1)
            temp=(1+Data_intest(k,:))'.*PortTemp.PRP_L1{k};
            PortTemp.PRP_L1{k+1}=temp./sum(temp);
        end
        re_PRP_L1=[];
        for k=1:size(Data_intest,1)
            re_PRP_L1=[re_PRP_L1 Data_intest(k,:)*PortTemp.PRP_L1{k}];  %*************************** PortNew1 should be replaced
        end
        return_PRP_L1 = sum(re_PRP_L1)/length(re_PRP_L1);
        vol_PRP_L1 = norm(re_PRP_L1-return_PRP_L1)/sqrt(length(re_PRP_L1)-1);
        %             SR_PRP_L1=[SR_PRP_L1 return_PRP_L1/vol_PRP_L1];
        %             SR_PRP_L1=[SR_PRP_L1 nu*return_PRP_L1-(vol_PRP_L1)^2];
        %             if length(Wealth.PRP_L1)==0
        %                 W0=1;
        %             else
        %                 W0=Wealth.PRP_L1(end);
        %             end
        money=cumprod(1+re_PRP_L1);
        SR_PRP_L1=[SR_PRP_L1 0.05*money(end)+nu*return_PRP_L1-(vol_PRP_L1)^2];
        if isempty(Lam1sL_best)
            Lam1sL_best=TOParam.lam1;
            Alpha2sL_best=TOParam.lam2;
        else if SR_PRP_L1(end)>SR_PRP_L1(end-1)
                Lam1sL_best=TOParam.lam1;
                Alpha2sL_best=TOParam.lam2;
            end
        end
    end
end

%% design the portfolio 
R=sqrtm(CovMatrix);
r_index=nu/2*(R\(MeanVec'));
TOParam.lam1=Lam1sL_best;
TOParam.lam2=Alpha2sL_best;
ScaleCovMatrix = 1e4 .* CovMatrix;
M = cell(NumStocks,1);
for whichi = 1:NumStocks
    tmpM = sparse(NumStocks,NumStocks);
    tmpM(whichi,:) = ScaleCovMatrix(whichi,:);
    M{whichi} = (tmpM + tmpM');
end
tau = 0.01*trace(ScaleCovMatrix)/(2*NumStocks);
w0 = oneN ./ NumStocks;
tmp = RP_Regularized_Linear(M, w0, tau, TOParam, ScaleCovMatrix, ApproxFun, 200, R, r_index);
w_PRP_L1 = tmp.w;