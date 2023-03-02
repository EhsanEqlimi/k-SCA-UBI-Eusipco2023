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
m_Set=[3]; % The number of sensors a.k.a the ambient dimension
n_Set=[10];% The number of sources (ith row correspond ith elemnet in m_se
% SNR_Add=[30]; % SNR values of additive noise
Sigma_Off=[1e-3]; % Parameter controls the standard deviation of normal noise over the active lements
T= 2000;%c*30; %1000; % The Number of data poins a.k.a the number of samples (time points)
AMode=1; % k-SCA Condition for A is satisfied
Iter_Num_A=1000;% The number of iteration to generate a good mixing matrix
Rank_Th_A=0.1; % to generate a good mixing matrix
Th_RansacSubspace=1e-4;% RANSAC Threshold in Distance Function for subspace estimation
ThBBC=1e-4; %BBC clustering for mixing idetification
Th_RansacMixing=1e-4;% RANSAC Threshold in Distance Function for mixing matrix identifiation
DualMode=0; %Twice subpspace clusteing
BAS_Set=[];NMSE_Set=[];Norm2Error_Set=[];Sub_Error_Set=[];Error_Set=[];
%% SCA Mixing (Simulted Sources)
%%'kSCA';%'kSCA';% Uniform Sparse Component Mixing  %'MSCA'=>Multiple Sparse Component Mixing, 'PermKSCA'
MixingMode='PermkSCA';% {'kSCANoisy'%'PermkSCA'%kSCA';%'MSCA';%OnOffGauss};
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
                X=X./(ones(m,1)*sqrt(R));
                X=FnColNormalizer(X);
                % Ahat=PKDSC(X, n, k, Sigma_B, Sigma_A, N_B, q, L_B, L_A, TH1, TH2, TH3);
                TrureConnMat=FnCreateConnMat(Labels);
                % figure,imagesc(TrureConnMat);title('True Connectivity Matrix');
                % X=FnColNormalizer(X);% It is impotant when wh have clos orgin samples in X;
                %% Subspace Estimation
                tic
                % % % % % %                 X=M;k=3;n=size(D1,2);c=nchoosek(n,k);
%                 [QRSubspaceInds,Clusters,SubSpaces,ComplementOrthofSubSpaces,ConnMat]=FnSubSpaceFindV3(X,Th_RansacSubspace,k,c,SubspaceInds,MaxIterEst,p);
%                 %% Mixing Matrix Identification.
%                 % Option 1
%                 % [Ahat,MinEV]=FnMixingCalc_Ehsan(ComplementOrthofSubSpaces,k,c,n); % Slow/High comp. cost/ needs all comb i.e. Order=Cx=nchoosek(1:c,g)
%                 % [Ahat,MinEV]=FnMixingCalc_EhsanNew(ComplementOrthofSubSpaces,k,c,n);
%                 %% Option 2 :RANSCA-based method
%                 [ComplementOrthofSubSpaces_A,Ahat]=FnSubSpaceFind_Mixing_Final(squeeze(ComplementOrthofSubSpaces),ThBBC,k,n,Th_RansacMixing,p,MaxIterEst);
%                 Time_Ours(Iter_Loop,n_Loop)=toc;
%                 %% Error estimation
%                 [ClusterinError1(Iter_Loop,n_Loop),ClusterinError2(Iter_Loop,n_Loop)]=FnSubspaceClusteringErrorFinder(Clusters,Labels);
%                 [BAS_Ours(Iter_Loop,n_Loop), MixingVectorerror_Ours{Iter_Loop,n_Loop},NMSE_Ours{Iter_Loop,n_Loop},...
%                     NMSSum_Ours(Iter_Loop,n_Loop),AhatNew_Ours{Iter_Loop,n_Loop},Norm2Error_Ours(Iter_Loop,n_Loop),Error_Ours(Iter_Loop,n_Loop),BAS_Detected_Ours(Iter_Loop,n_Loop)]...
%                     = FnMixingIdentificationError(A,Ahat);
%                 disp(['BAS_Ours=' num2str(BAS_Ours(Iter_Loop,n_Loop))])
%                 disp(['BAS_Detected_Ours=' num2str(BAS_Detected_Ours(Iter_Loop,n_Loop))])
                %%
                tic
                AhatMovahedi=FnMovahediUBI(m,n,k,T,X);
                Time_Movahedi(Iter_Loop,n_Loop)=toc;
                size(AhatMovahedi)
                [BAS_Movahedi(Iter_Loop,n_Loop), MixingVectorerror_Movahedi{Iter_Loop,n_Loop},NMSE_Movahedi{Iter_Loop,n_Loop},NMSSum_Movahedi(Iter_Loop,n_Loop),...
                    AhatNew_Movahedi{Iter_Loop,n_Loop},Norm2Error_Movahedi(Iter_Loop,n_Loop),Error_Movahedi(Iter_Loop,n_Loop),BAS_Detected_Movahedi(Iter_Loop,n_Loop)] = ...
                    FnMixingIdentificationError(A,AhatMovahedi);
                A_Set{Iter_Loop,n_Loop}=A;
                S_Set{Iter_Loop,n_Loop}=S;
                X_Set{Iter_Loop,n_Loop}=X;
%                 OCS_Set{Iter_Loop,n_Loop}=ComplementOrthofSubSpaces;
%                 disp(['BAS_Movahedi=' num2str(BAS_Movahedi(Iter_Loop,n_Loop))])
%                 disp(['BAS_Detected_Movahedi=' num2str(BAS_Detected_Movahedi(Iter_Loop,n_Loop))])
            end
        end
    end
end
% % end
%  load  Ransa_Movahadi_100trial.mat
% for i=1:size(MixingVectorerror_Movahedi,1)
%     for j=1:1:size(MixingVectorerror_Movahedi,2)
%         nhat_mov(i,j)=length(find(MixingVectorerror_Movahedi{i,j}<.1));
% %           nhat_ours(i,j)=length(find(MixingVectorerror_Ours{i,j}<.1));
% 
%     
%     end
% end
% % 
% % b = bar(1:10,[mean(nhat_mov,1) mean(nhat_ours,1)]);
% % b.FaceColor = 'g';
% % figure,stem(1:5,mean(nhat_mov,1),'--k>')
% % hold on
% % stem(1:5,mean(nhat_ours,1),'-*rs')
% % bpcombined = [mean(nhat_mov,1), mean(nhat_ours,1)];
% % hb = bar(1:10, bpcombined, 'grouped');
% close all
% XTick={'[3,5]', '[3,6]','[3,8]','[3,9]','[3,10]'};
% 
% % title('Blocking Probability vs Routing Level');
% Legend={'\sigma_{off} = 10^{-6}','\sigma_{off}=10^{-4}','\sigma_{off}=10^{-3}'};
% % plot(BAS)
% %% Plot Config.
% 
% figure,
% plot(1:5,log10(mean(BAS_Detected_Ours,1)),'-.rs','MarkerSize',8,'LineWidth',2);
% LogMeanBAS_Ours=log10(mean(BAS_Detected_Ours,1));
% LogMeanBAS_Mov=log10(mean(BAS_Detected_Movahedi,1));
% 
% xlim([.85 5.32])
% 
% text([1 1.6 2.7 3.7 4.7],[-2.9,-2.78, -2.4, -2.24, -2.10 ] , '$$\hat{n}$$','Interpreter','Latex','fontsize', 17,'fontweight','bold' )
% % text([1.1 1.8 3 4 4.9],LogMeanBAS_Ours+.1 , '$$\hat{n}$$','Interpreter','Latex','fontsize', 17,'fontweight','bold' )
% 
% hold on
% plot(1:5,log10(mean(BAS_Detected_Movahedi,1)),'--k>','MarkerSize',8,'LineWidth',2);
% text([1.1 1.8 3 4 4.9],LogMeanBAS_Mov+.1 , '$$\hat{n}$$','Interpreter','Latex','fontsize', 16,'fontweight','bold' )
% 
% set(gca, 'XTick', 1:5,...
%     'XTickLabel', XTick);
% % xticklabels(XTick,'Interpreter','Latex' )
% % set(gca, 'fontsize', 10,'fontweight','bold','MarkerSize',15);
% 
% %    set(gca,'YTick',[-pi 0 pi], 'YTickLabel', {'-\pi','0','\pi'}, 'fontsize', 18);
% %     set(gca, 'XTickLabel', , 'fontsize', 6);
% 
% xlabel('$$[m , n]$$','Interpreter','Latex','fontsize', 16,'fontweight','bold' );
% ylabel('$$\log_{10}(\textrm{BAS})$$','Interpreter','Latex','fontsize', 16,'fontweight','bold' )
% legend({'Proposed Algorithm','Algorithim in [16]'},'location','Best')
% set(gca, 'fontsize', 10,'fontweight','bold');
% box on
% 
% figure,
% data = [mean(nhat_ours,1)' , mean(nhat_mov,1)' ];
% hb = bar(data);
% set(hb(1), 'FaceColor','r')
% set(hb(2), 'FaceColor','k')
% ylabel('Average number of identified vectors');
% xlabel('[m,n]');
% legend({'Proposed Algorithm','Algorithim in [16]'},'location','Best')
% set(gca, 'XTick', 1:5,...
%     'XTickLabel', XTick,'Interpreter','Latex','fontsize', 16,'fontweight','bold' );
% set(gca, 'fontsize', 10,'fontweight','bold');
% box on
% %%
% figure,open( 'Fig2_text.fig')
% xlabel('$$[m , n]$$','Interpreter','Latex','fontsize', 16,'fontweight','bold' );
% ylabel('$$\log_{10}(\textrm{BAS})$$','Interpreter','Latex','fontsize', 16,'fontweight','bold' )