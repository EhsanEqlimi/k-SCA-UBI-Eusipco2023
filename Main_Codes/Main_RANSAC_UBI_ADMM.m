%**************************************************************************
% Main file to run the EEG microstates based on k-Sparse Componenet
% Analysis
%Ehsan Eqlimi Oct 2020
%**************************************************************************
clc;
clear ;
close all;
warning off

%% Preliminaries
Iter_Num=100; % The number of iteration
m_Set=[3,4,5]; % The number of sensors a.k.a the ambient dimension
n_Set=[4,5,6,8;5,6,7,8;6,7,8,9];% The number of sources (ith row correspond ith elemnet in m_set
m_Set=[64]; % The number of sensors a.k.a the ambient dimension
n_Set=[4];% The number of sources (ith row correspond ith elemnet in m_set

% SNR_Add=[30]; % SNR values of additive noise
Sigma_Off=[0]; % Parameter controls the standard deviation of normal noise over the active lements

T= 1000;%c*30; %1000; % The Number of data poins a.k.a the number of samples (time points)
% zero sources (non-strict saprsity)
AMode=1; % k-SCA Condition for A is satisfied
Iter_Num_A=1000;% The number of iteration to generate a good mixing matrix
Rank_Th_A=0.1; % to generate a good mixing matrix
Th_RansacSubspace=1e-5;% RANSAC Threshold in Distance Function for subspace estimation
ThBBC=1e-3; %BBC clustering for mixing idetification
Th_RansacMixing=1e-5;% RANSAC Threshold in Distance Function for mixing matrix identifiation

DualMode=0; %Twice subpspace clusteing
BAS_Set=[];NMSE_Set=[];Norm2Error_Set=[];Sub_Error_Set=[];Error_Set=[];
%% SCA Mixing (Simulted Sources)
%%'kSCA';%'kSCA';% Uniform Sparse Component Mixing  %'MSCA'=>Multiple Sparse Component Mixing, 'PermKSCA'
MixingMode='kSCA';% {'kSCANoisy'%'PermkSCA'%kSCA';%'MSCA';%OnOffGauss};
Orth=0; %if Orth=1 A is orth
p=0.95;
% Sigma=0;
CT=0;
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
            k=2;
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
                %                 g=sqrt(sum(X.^2,1));
                %                 X=X(:,g>.2);
                %                 Labels=Labels(g>.2);
                %                 %normalizing every column of the mixture matrix
                %                 R=sum(X.^2,1);
                %                 X=X./(ones(m,1)*sqrt(R));
                %                 X=FnColNormalizer(X);
                % Ahat=PKDSC(X, n, k, Sigma_B, Sigma_A, N_B, q, L_B, L_A, TH1, TH2, TH3);
                TrureConnMat=FnCreateConnMat(Labels);
                % figure,imagesc(TrureConnMat);title('True Connectivity Matrix');
                % X=FnColNormalizer(X);% It is impotant when wh have clos orgin samples in X;
                ColorLab={'r','b','g','m'};
                
                figure,
                for i=1:n
                    %     Temp=[zeros(1,(i-1)*L*OL) Cosine zeros(1,n*L-L-(i-1)*L*OL)] ;%figure,plot(S);
                    %     S(i,:)=Temp(:,1:RealL);
                    subplot(n,1,i)
                    plot(1:size(S,2),S(i,:),ColorLab{i})
                end
                %% Subspace Estimation
                tic
                % % % % % %                 X=M;k=3;n=size(D1,2);c=nchoosek(n,k);
                Xhat=X-repmat(mean(X,1),size(X,1),1);
                [Coeff,Score,Latent,TSquared] =pca((X'));
                % Data Projection based on PCA
                %                 Xp = FnDataProjection(X,r,'PCA');
                Colors={'.r', '.b', '.g', '.m','-y','..r'}; % Colors for clusters
                
                Actualc=4;
                [LabelsKmenas,KmeansClusters]=kmeans(X',Actualc);
                [~,StartIndClusters,~]=unique(LabelsKmenas, 'stable');
                StartIndClusters(4)=size(S,2);
                %                 figure,
                %                 for j=1:n
                %                     scatter3((StartIndClusters(j):StartIndClusters(j+1)-1,1),Score(StartIndClusters(j):StartIndClusters(j+1)-1,2),...
                %                         Score(StartIndClusters(j):StartIndClusters(j+1)-1,3),Colors{j});
                %                     hold on
                %                 end
                [QRSubspaceInds,Clusters,SubSpaces,ComplementOrthofSubSpaces,ConnMat]=FnSubSpaceFindV3(X,Th_RansacSubspace,k,c,SubspaceInds,MaxIterEst,p);
                %                 [inliers,M]=FnDistanceClacBetweenSubspaceVector(ComplementOrthofSubSpaces(:,:,1),X(:,1),.1)
                for SubInd=1:size(ComplementOrthofSubSpaces,3)
                    for DataInd=1:size(X,2)
                        Dist(SubInd,DataInd)=norm(ComplementOrthofSubSpaces(:,:,SubInd)'*X(:,DataInd));
                    end
                end
                %% Mixing Matrix Identification.
                % Option 1
                [Ahat,MinEV]=FnMixingCalc_Ehsan(ComplementOrthofSubSpaces,k,c,n); % Slow/High comp. cost/ needs all comb i.e. Order=Cx=nchoosek(1:c,g)
                %                                 [Ahat,MinEV]=FnMixingCalc_EhsanNew(ComplementOrthofSubSpaces,k,c,n);
                %% Option 2 :RANSCA-based method
                %                                 [ComplementOrthofSubSpaces_A,Ahat]=FnSubSpaceFind_Mixing_Final(squeeze(ComplementOrthofSubSpaces),ThBBC,k,n,Th_RansacMixing,p,MaxIterEst);
                ET(CT)=toc;
                %% Error estimation
                [ClusterinError1,ClusterinError2]=FnSubspaceClusteringErrorFinder(Clusters,Labels);
                figure,stem(Clusters);title('my K-SCA');
                figure,stem(Labels);title('Ground truth');
                Missrate = Misclassification(Clusters',Labels')
                
                
                
                [MixingIdentificatioError, MixingVectorerror,NMSE,NMSSum,AhatNew,Norm2Error,Error] = FnMixingIdentificationError(A,Ahat);
                
                [LabelsKmenas,KmeansClusters]=kmeans(Xhat',n);
                [MixingIdentificatioError2, MixingVectorerror2,NMSE2,NMSSum2,AhatNew2,Norm2Error2,Error2] = FnMixingIdentificationError(A,KmeansClusters');
                
                load TempalteEEG.mat;
                EEG.msinfo.MSMaps(4).Maps=AhatNew2';
                pop_ShowIndMSMaps(EEG,4);
                
                EEG.msinfo.MSMaps(4).Maps=AhatNew';
                pop_ShowIndMSMaps(EEG,4);
                
                
                EEG.msinfo.MSMaps(4).Maps=A';
                pop_ShowIndMSMaps(EEG,4);
                %                                 ChanLocs = readlocs('ChanLocsEE.loc'); % Read 64 EEG channels (Biosemi 10-10);
                %
                %                                 for mm=1:4
                %                                     figure,
                %                                     topoplot(AhatNew(:,mm), ChanLocs,'electrodes','on','style','map','plotrad',.54);
                %                                 end
                %
                %                                 for mm=1:4
                %                                     figure,
                %                                     topoplot(A(:,mm), ChanLocs,'electrodes','on','style','map','plotrad',.54);
                %                                 end
                %                                 for mm=1:4
                %                                     figure,
                %                                     topoplot(AhatNew2(:,mm), ChanLocs,'electrodes','on','style','map','plotrad',.54);
                %                                 end
                
                
                %%
                %                 tic
                %                 AhatMovahedi=FnMovahediUBI(m,n,k,T,X);
                %                 Time(Iter_Loop)=toc;
                %
                %                 [MixingIdentificatioError, MixingVectorerror{Iter_Loop},NMSE{Iter_Loop},NMSSum,AhatNew,Norm2Error,Error] = FnMixingIdentificationError(A,AhatMovahedi);
                %
                % %                 Sub_Error_Set=[Sub_Error_Set ClusterinError1];
                %                 BAS_Set=[BAS_Set MixingIdentificatioError]
                %                 NMSE_Set=[NMSE_Set NMSSum];
                %                 Norm2Error_Set=[Norm2Error_Set Norm2Error];
                %                 Error_Set=[Error_Set Error];
                %
                %                 A_Set{CT}=A;
                %                 Ahat_Set{CT}=AhatMovahedi;
                %                 S_Set{CT}=S;
                %                 X_Set{CT}=X;
                %                 AhatNew_Set{CT}=AhatNew;
                %% ADMM-SSC
                Actualc=c;
                alpha=60;
                r = 0; affine = false; outlier = true; rho = 1;
                [missrate,C,grps,Clus] = SSC(X,r,affine,alpha,outlier,rho,Labels);
                Sub=FnSubspaceCalc(X,grps',1,Actualc);
                Subb=squeeze(Sub)
                figure,stem(grps)
                title('SCA-ADMM-Alpha-10')
                
                [Sub,OrthSub]=FnSubspaceCalc(X,grps',2,6);
                Subb=squeeze(Sub);
                % Ahat=FnMixingCalc(OrthSub,2,6,4)
                [Ahat,MinEV]=FnMixingCalc_Ehsan(OrthSub,2,4,4)
                [MixingIdentificatioError, MixingVectorerror,NMSE,NMSSum,AhatNew,Norm2Error,Error] = FnMixingIdentificationError(A,Ahat)
                
                
                 for SubInd=1:size(OrthSub,3)
                    for DataInd=1:size(X,2)
                        DistADMM(SubInd,DataInd)=norm(OrthSub(:,:,SubInd)'*X(:,DataInd));
                    end
                end
                figure,stem(DistADMM')
                
                %%K-means
                %%K-means
                [LabelsKmenas,KmeansClusters]=kmeans(Xhat',Actualc);
                figure,stem(LabelsKmenas)
                title('k-means')
                Missrate_Kmeans = Misclassification(LabelsKmenas,GroundTruthLables)
            end
        end
    end
end
% end


