%% clear the environment
clear all; clc;

%% load the data & define new varaibles
% load the data
load('RETURN_AllAssets_SP100_20161025.mat');
RETURN(:,1)=[];
NumStocks=size(RETURN,2);

% define the number of paris of in-sample and uot-of-sample windows
NumWindows=27;  % set by the designer
for i=1:NumWindows
    RollingData(i).return=RETURN(60*(i-1)+1:183+60*(i-1+1),:) ; 
end

% define the variable to store the results
PRP_L_Laqmbda=[];  % store the slected value of lambda1 for linear approxiamtion 
PRP_L_Alpha=[];   % store the slected value of lambda2 for linear approxiamtion 
PRP_Q_Laqmbda=[]; % store the slected value of lambda1 for quadratic approxiamtion 
PRP_Q_Alpha=[]; % store the slected value of lambda2 for quadratic approxiamtion 

nu=1; % the parameter for (risk-nu*mean) 

%% portfolio design of each window
for i=1:NumWindows
    %% divide rollingwindows
    % the rolling window = an in-sample window + an out-of-sample window. 
    % For the cross validation, the in-sample window = an in-sample training window + an in-sample testing window)
    [Data_p, Data_in, Data_out, Data_intrain, Data_intest, Data_out_p]=DivideWindows(RollingData(i).return);  
    %% estimate the covariance matrices and mean vectors
    [CovMatrix, CovMatrix_train, CovMatrix_test, CovMatrix1,MeanVec, MeanVec_train]=EstimateCovarianceMean(Data_in,Data_intrain,Data_intest,Data_out_p,Data_p); 
    
     %% Sandard Mean-varaince portfolio 
%     % use MOSEK
%     [Port(i).MV]=MVPort(CovMatrix,MeanVec,nu);   
%     %% Sparse MV portfolio (L2+L0 regularization: Greedy algorithm)
%     NumSlc=[4 7 7 3 5 4 2 5 4 4 4 8 7 9 7 6 8 25 29 13 14 15 10 9 4 5 8];
%     % use MOSEK
%     tau=0.01;
%     [Port(i).LG]=SparseMVPort(CovMatrix,MeanVec,tau,nu,NumSlc(i)) ;
%     %% EW portfolio
%     oneN = ones(NumStocks, 1);
%     Port(i).EW = oneN ./ NumStocks;  
%     %% ERC portfolio
%     MaxIter=100;
%     [Port(i).ERC]=ERCPort(oneN,NumStocks,CovMatrix,MaxIter);  

    %% Proposed portfolio with quadratic approximation
    Lam1s=[0.00001:0.00001: 0.00001];
    Lam2s=[0.001:0.001:0.001];
    [Port(i).PRP, Lam1sQ_best, Alpha2sQ_best]=RiskParityPort_Quadratic(Data_intest,CovMatrix_train,MeanVec_train,CovMatrix,MeanVec,Lam1s,Lam2s,nu,NumStocks);
    PRP_Q_Laqmbda=[PRP_Q_Laqmbda Lam1sQ_best];
    PRP_Q_Alpha=[PRP_Q_Alpha Alpha2sQ_best];

    %% Proposed portfolio with linear approximation
    Lam1s=[0.00001:0.00001: 0.00001];
    Lam2s=[0.001:0.001:0.001];
    [Port(i).PRP_L1, Lam1sL_best, Alpha2sL_best]=RiskParityPort_Linear(Data_intest,CovMatrix_train,MeanVec_train,CovMatrix,MeanVec,Lam1s,Lam2s,nu,NumStocks);
    PRP_L_Laqmbda=[PRP_L_Laqmbda Lam1sL_best];
    PRP_L_Alpha=[PRP_L_Alpha Alpha2sL_best];
    
    %% Caculate the portfolio weights for each trading day of the out-of-sample window
    % since we don't reblance everyday, the portfolio vector changes everyday. We record all these vectors for further analysis.
    
    % Portfolio weights in the out-of-sample window (exclude the weekends)
    PortNew1(i).EW{1}=Port(i).EW;
    PortNew1(i).ERC{1}=Port(i).ERC;
    PortNew1(i).MV{1}=Port(i).MV;
    PortNew1(i).LG{1}=Port(i).LG;
    PortNew1(i).PRP{1}=Port(i).PRP;
    PortNew1(i).PRP_L1{1}=Port(i).PRP_L1;
    TEMP=1+Data_out_p;
    for k=1:size(TEMP,1)
        temp=TEMP(k,:)'.*PortNew1(i).EW{k};
        PortNew1(i).EW{k+1}=temp./sum(temp);
        temp=TEMP(k,:)'.*PortNew1(i).ERC{k};
        PortNew1(i).ERC{k+1}=temp./sum(temp);
        temp=TEMP(k,:)'.*PortNew1(i).MV{k};
        PortNew1(i).MV{k+1}=temp./sum(temp);
        temp=TEMP(k,:)'.*PortNew1(i).LG{k};
        PortNew1(i).LG{k+1}=temp./sum(temp);
        temp=TEMP(k,:)'.*PortNew1(i).PRP{k};
        PortNew1(i).PRP{k+1}=temp./sum(temp);
        temp=TEMP(k,:)'.*PortNew1(i).PRP_L1{k};
        PortNew1(i).PRP_L1{k+1}=temp./sum(temp);
    end
    
    % Portfolio weights in the out-of-sample window (include the weekends)
    PortNew(i).EW{1}=Port(i).EW;
    PortNew(i).ERC{1}=Port(i).ERC;
    PortNew(i).MV{1}=Port(i).MV;
    PortNew(i).LG{1}=Port(i).LG;
    PortNew(i).PRP{1}=Port(i).PRP;
    PortNew(i).PRP_L1{1}=Port(i).PRP_L1;
    TEMP=1+Data_out;
    for k=1:size(TEMP,1)
        temp=TEMP(k,:)'.*PortNew(i).EW{k};
        PortNew(i).EW{k+1}=temp./sum(temp);
        temp=TEMP(k,:)'.*PortNew(i).ERC{k};
        PortNew(i).ERC{k+1}=temp./sum(temp);
        temp=TEMP(k,:)'.*PortNew(i).MV{k};
        PortNew(i).MV{k+1}=temp./sum(temp);
        temp=TEMP(k,:)'.*PortNew(i).LG{k};
        PortNew(i).LG{k+1}=temp./sum(temp);
        temp=TEMP(k,:)'.*PortNew(i).PRP{k};
        PortNew(i).PRP{k+1}=temp./sum(temp);
        temp=TEMP(k,:)'.*PortNew(i).PRP_L1{k};
        PortNew(i).PRP_L1{k+1}=temp./sum(temp);
    end
end

%% save the results
save('DesignedPorts.mat','NumWindows','nu','RollingData','Port','PortNew1','PortNew');
