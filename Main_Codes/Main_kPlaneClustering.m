%**************************************************************************
% Main and demo m-file to run on-line k-plane clustering
%**************************************************************************
%**************************************************************************
% Main file to run Gram-Shcmidt and RANSAC based k-SCA UBI
% 20/12/94
% Ehsan Eqlimi,TUMS
%Last edit: 3/02/2020 by EE
%**************************************************************************
clc;
clear ;
close all;
warning off

%% Preliminaries
Iter_Num=100; % The number of iteration
m_Set=[3,4,5]; % The number of sensors a.k.a the ambient dimension
n_Set=[4,5,6,8;5,6,7,8;6,7,8,9];% The number of sources (ith row correspond ith elemnet in m_set
m_Set=[3]; % The number of sensors a.k.a the ambient dimension
n_Set=[5];% The number of sources (ith row correspond ith elemnet in m_set

% SNR_Add=[30]; % SNR values of additive noise
Sigma_Off= 0;%[1e-4]; % Parameter controls the standard deviation of normal noise over the active lements

T= 2000;%c*30; %1000; % The Number of data poins a.k.a the number of samples (time points)
% zero sources (non-strict saprsity)
AMode=1; % k-SCA Condition for A is satisfied
Iter_Num_A=1000;% The number of iteration to generate a good mixing matrix
Rank_Th_A=0.1; % to generate a good mixing matrix
Th_RansacSubspace=1e-5;% RANSAC Threshold in Distance Function for subspace estimation
ThBBC=1e-4; %BBC clustering for mixing idetification
Th_RansacMixing=1e-5;% RANSAC Threshold in Distance Function for mixing matrix identifiation

DualMode=0; %Twice subpspace clusteing
BAS_Set=[];NMSE_Set=[];Norm2Error_Set=[];Sub_Error_Set=[];Error_Set=[];
%% SCA Mixing (Simulted Sources)
%%'kSCA';%'kSCA';% Uniform Sparse Component Mixing  %'MSCA'=>Multiple Sparse Component Mixing, 'PermKSCA'
MixingMode='PermkSCA';% {'kSCANoisy'%'PermkSCA'%kSCA';%'MSCA';%OnOffGauss};
Orth=0; %if Orth=1 A is orth
p=0.95;
% Sigma=0;
CT=0;
LearningTime=20;
%% Main Process
for Iter_Loop=1:Iter_Num
    disp(Iter_Loop);
    %     for SNR_Add_Loop=1:length(SNR_Add)
    %         SNRIn=SNR_Add(SNR_Add_Loop);
    for Sigma_Off_Loop=1:length(Sigma_Off)
        Sigma=Sigma_Off(Sigma_Off_Loop);
        
        for m_Loop=1:length(m_Set)
            CT=CT+1;
            m=m_Set(m_Loop);
            k=m-1; % The number of acitive source in each time point (nonzero elements
            %in each col of S) or the dimension of each subspace
            r=m-k; % The dimension of the normal space (normal vector/orthogonal complement) of each subspace
            n_Correpond=n_Set(m_Loop,:);
            
            for n_Loop=1:size(n_Set,2)
                n=n_Correpond(n_Loop);
                c=nchoosek(n,k); % The number of all possible r-dim subspaces
                N=zeros(1,c);
                N(1,:)=ceil(T/c); % N is the number of the subspaces in k-SCA mode
                for j=1:k
                    for i=1:nchoosek(n,j)
                        Nk(i,j)=ceil(ceil(T/k)/nchoosek(n,j));% N showes the number of each subspace in MSCA mode
                    end
                end
                e= 1-(1/c);%(c-1)/c;  % Probability that a point is an outlier
                %                     MaxIterEst=log(1-p)/log(1-(1-e)^k); %1e6;%log(1-p)/log(1-(1-e)^k);
                MaxIterEst=1e6;
                
                if MaxIterEst<1000% Sometimes if k=1, it happnes
                    MaxIterEst=1e4;
                end
                %% Sparse Component Mixing (Simulation)
                [X,S,OrthA,A,Labels,SubspaceInds,SubspaceNum]=FnSparseComponentMixing(m,n,k,T,N,Nk,Sigma,Iter_Num_A,Rank_Th_A,AMode,MixingMode);
                % Add noise
                %                     for i=1:size(X,1)
                %                         SumSquareX=sum(X(i,:).^2);
                %                         [SigmaIn(i),VarIn(i)]=FnSNR2Sigma(SumSquareX,SNRIn,size(X,2));
                %                         Noise(i,:)=SigmaIn(i).*randn(1,size(X,2));
                %                         XNoisy(i,:)=X(i,:)+ Noise(i,:);
                %                         SNRIn_2(i)=10*log10(SumSquareX/(sum( Noise(i,:).^2)));
                %                     end
                %                     Th_RansacSubspace1=Sigma ;% RANSAC Threshold in Distance Function for subspace estimation
                %                     Th_RansacSubspace2=mean(SigmaIn) ;% RANSAC Threshold in Distance Function for subspace estimation
                %                     Th_RansacSubspace=0.1*max(Th_RansacSubspace1,Th_RansacSubspace2);
                %                 Th_RansacSubspace=1e-6;
                % %
                %
                %                     X=XNoisy;
                %Omitting the points that are near the origin
                g=sqrt(sum(X.^2,1));
                X=X(:,g>.2);
                Labels=Labels(g>.2);
                %normalizing every column of the mixture matrix
                R=sum(X.^2,1);
                %                 X=X./(ones(m,1)*sqrt(R));
                %                 X=FnColNormalizer(X);
                % Ahat=PKDSC(X, n, k, Sigma_B, Sigma_A, N_B, q, L_B, L_A, TH1, TH2, TH3);
                TrureConnMat=FnCreateConnMat(Labels);
                % figure,imagesc(TrureConnMat);title('True Connectivity Matrix');
                X=FnColNormalizer(X);% It is impotant when wh have clos orgin samples in X;
                %                 %% Washizawa
                %                 Ct=0;
                %                 SelectedIndices=1:size(X,2);
                %                 for j=1:c
                %                     disp(j)
                %                     W(:,j)=FnColNormalizer(rand(m,1));
                %                     OutliersInds=[];Ct=0;
                %                     for t=1:LearningTime*95
                %                         %                         disp(t);
                %                         for i=SelectedIndices
                %                             Eta(t)=(0.1- (5*1e-5)*t);
                %                             Theta(t)=cos(pi/4+(4e-4)*t);
                %                             if(abs(dot(W(:,j),X(:,i))))<=Theta(t) %theta(t)
                %                                 W(:,j)=W(:,j)-Eta(t).*dot(W(:,j),X(:,i)).*X(:,i);
                %                                 %                                 W(:,j)=FnColNormalizer(W(:,j));
                %                                 W(:,j)= W(:,j)/norm( W(:,j));
                %                             end
                %
                %                         end
                %                     end
                %                     % remove fitting samples
                %                     distance = X(:,SelectedIndices)' * W(:,j);
                %                     Outliers=find((abs(distance)) <= Theta(t));
                %                     SelectedIndices=setdiff(SelectedIndices,(Outliers'));
                
                
                
                
                
                
                
                %                     OutliersInds=[];
                %                     Ct=0;
                %                     for i=SelectedIndices
                %                         if abs(dot(W(:,j),X(:,i)))<= Theta(t) %theta(t)
                %                             Ct=Ct+1;
                %                             OutliersInds(Ct)=i;
                %
                %                         end
                %                     end
                
                
                %                     SelectedIndices=OutliersInds;
                
            end
            %                 W
            %                 [QRSubspaceInds,Clusters,SubSpaces,ComplementOrthofSubSpaces,ConnMat]=FnSubSpaceFindV3(X,Th_RansacSubspace,k,c,SubspaceInds,MaxIterEst,p);
            %                 ComplementOrthofSubSpaces
            %                 %                 learntimes=500;norm_th=0;th_init=cos(pi/4+(4e-4)*1);th_dec=1e-5;th_min=cos(pi/4+(4e-4)*100);eta_ini=(0.1- (5*1e-5)*1);eta_dec=0.1- (5*1e-5);eta_min=(0.1- (5*1e-5)*100);
            %                 normalvectors=GetNormalVectors(X,learntimes,norm_th,th_init,th_dec,th_min,eta_ini,eta_dec,eta_min,c);
            %                 normalvectors
            
            %                 Ahat=FnMixingCalc_Ehsan(W,k,c,n);
            A_Set=
            tic
            Ahatmov=FnMovahediUBI(m,n,k,T,X);
            Time(Iter_Loop)=toc;
            Time(Iter_Loop)
            [MixingIdentificatioError(Iter_Loop), MixingVectorerror{Iter_Loop},NMSE{Iter_Loop},NMSSum{Iter_Loop},AhatNew{Iter_Loop},Norm2Error(Iter_Loop),Error(Iter_Loop)] =FnMixingIdentificationError(A,Ahatmov)
        end
    end
end
% end
