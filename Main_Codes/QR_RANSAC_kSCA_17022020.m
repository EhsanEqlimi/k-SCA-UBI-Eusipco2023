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
m=6; % The number of sensors a.k.a the ambient dimension
n=12; % The number of sources
k=2; % The number of acitive source in each time point (nonzero elements
%in each col of S) or the dimension of each subspace
r=m-k; % The dimension of the normal space (normal vector/orthogonal complement) of each subspace
c=nchoosek(n,k); % The number of all possible r-dim subspaces
T= c*m*10; %1000; % The Number of data poins a.k.a the number of samples (time points)
Sigma=1e-3; % Parameter controls the standard deviation of normal noise over the
% zero sources (non-strict saprsity)
AMode=1; % k-SCA Condition for A are satisfied
IterNum=500;% The number of iteration to generate a good mixing matrix
RankTh=0.1; % to generate a good mixing matrix
Th_RansacSubspace=1e-4;% RANSAC Threshold in Distance Function for subspace estimation
Th3=1e-8;
Sigma1=0.02;
Th_RansacMixing=1e-5;% RANSAC Threshold in Distance Function for mixing matrix identifiation
Iter_MixingId=50*c; % is not importanat!
p=0.95; % Desired probability that we get a good sample in ransac
e= 1-(1/c);%(c-1)/c;  % Probability that a point is an outlier
MaxIterEst= log(1-p)/log(1-(1-e)^k); %1e6;%log(1-p)/log(1-(1-e)^k);
if MaxIterEst<1000
MaxIterEst=1e5;
end
%% RANSAC Functions:
DistFunc = @FnDistanceClacBetweenSubspaceVector;
DegFunc= @(x) FnDegenerateEval(x,k);
FitFunc=@(x) FnPNVClac(x);

g=nchoosek(n-1,k-1);
DegFunc2= @(x) FnDegenerateEval(x,m-1);

%% SCA Mixing (Simulted Sources)
%%'kSCA';%'kSCA';% Uniform Sparse Component Mixing  %'MSCA'=>Multiple Sparse Component Mixing, 'PermKSCA'
MixingMode='PermkSCA';% {'kSCANoisy'%'PermkSCA'%kSCA';%'MSCA';%OnOffGauss};
Orth=0; %if Orth=1 A is orth
N=zeros(1,c);
N(1,:)=ceil(T/c); % N is the number of the subspaces in k-SCA mode
%NK: The cardinal number of  each subspace in MSCA mode (only used for
%MSCA)
Nk=zeros(c,k);
for j=1:k
    for i=1:nchoosek(n,j)
        Nk(i,j)=ceil(ceil(T/k)/nchoosek(n,j));% N showes the number of each subspace in MSCA mode
    end
end

%% Sparse Component Mixing (Simulation)
[X,S,OrthA,A,Labels,SubspaceInds,SubspaceNum]=FnSparseComponentMixing(m,n,k,T,N,Nk,Sigma,IterNum,RankTh,AMode,MixingMode);
TrureConnMat=FnCreateConnMat(Labels);
figure,imagesc(TrureConnMat);title('True Connectivity Matrix')

X=FnColNormalizer(X);
% figure, scatter(X(1,:),X(2,:));
% X=X+0.001*randn(m,T);
% [Sjj]=L1_norm_min(X,A);
%% Subspace Estimation
% % % % % % % % % [IDX SS] = ksubspaces(X',c,k);
% % % % % % % % % [aIDX aSS] = seqksubspaces(X,c,0.25);
% % %  CMat = admmOutlier_mat_func(X,1,20);
% % %  CKSym = BuildAdjacency(thrC(CMat,.1));
% % % grps = SpectralClustering(CKSym,c);
% % % grps = bestMap(Labels,grps');
% % % missrate = sum(Labels(:) ~= grps(:)) / length(Labels);
[QRSubspaceInds,Clusters,SubSpaces,ComplementOrthofSubSpaces,ConnMat]=FnSubSpaceFind(X,Th_RansacSubspace,k,c,SubspaceInds,DistFunc,DegFunc,FitFunc,MaxIterEst,p);
figure,imagesc(ConnMat);title('Estimated Connectivity Matrix')
[ClusterinError1,ClusterinError2]=FnSubspaceClusteringErrorFinder(Clusters,Labels)

% [QRSubspaceInds,Clusters,SubSpaces,ComplementOrthofSubSpaces,ConnMat]=FnSubSpaceFindV2(X,Th_RansacSubspace,k,c,SubspaceInds,DistFunc,DegFunc,FitFunc);
% [vv,mm]=kmeans(ComplementOrthofSubSpaces',35);
% Clusters2=Clusters;
% ConnMat2=zeros(size(ConnMat));
% for tt=1:length(vv)
% Clusters2(find(Clusters==tt))=vv(tt);
% ConnMat2(find(Clusters2==vv(tt)),find(Clusters2==vv(tt)))=1;
% 
% end
% figure,imagesc(ConnMat2)
% figure,imagesc(ConnMat)
% [ClusterinError1,ClusterinError2]=FnSubspaceClusteringErrorFinder(IDX',Labels)


% Inliers=find(Clusters);
% [QRSubspaceInds,Clusters,SubSpaces,ComplementOrthofSubSpaces,ConnMat]=FnSubSpaceFind(X(:,Inliers),Th_RansacSubspace,k,c,SubspaceInds,DistFunc,DegFunc,FitFunc);
% 
% ClusterinError=FnSubspaceClusteringErrorFinder(Clusters,Labels)
% figure,imagesc(ConnMat)
%% Mixing Matrix Identification.
%% Option 1 
% [Ahat,MinEV]=FnMixingCalc_Ehsan(ComplementOrthofSubSpaces,k,c,n); % Slow/High comp. cost/ needs all comb i.e. Order=Cx=nchoosek(1:c,g)
% [Ahat,MinEV]=FnMixingCalc_EhsanNew(ComplementOrthofSubSpaces,k,c,n); 
%% Option 2 :RANSCA-based method


% [QRSubspaceInds2,Clusters2,SubSpaces2,ComplementOrthofSubSpaces2,ConnMat2,Ahat]=FnSubSpaceFind_Mixing(ComplementOrthofSubSpaces,Th_RansacMixing,k,Iter_MixingId,SubspaceInds,DistFunc,DegFunc,FitFunc,n,A);
% [QRSubspaceInds2,Clusters2,SubSpaces2,ComplementOrthofSubSpaces2,ConnMat2,Ahat]=FnSubSpaceFind_Mixing3(ComplementOrthofSubSpaces,Th_RansacMixing,k,Iter_MixingId,SubspaceInds,DistFunc,DegFunc,FitFunc,n,A);
[QRSubspaceInds2,Clusters2,SubSpaces2,ComplementOrthofSubSpaces2,ConnMat2,Ahat]=FnSubSpaceFind_Mixing5(ComplementOrthofSubSpaces,Th_RansacMixing,k,Iter_MixingId,SubspaceInds,DistFunc,DegFunc2,FitFunc,n,A);
% [QRSubspaceInds2,Clusters2,SubSpaces2,ComplementOrthofSubSpaces2,ConnMat2,Ahat]=FnSubSpaceFind_MixingNew(ComplementOrthofSubSpaces,Th_RansacMixing,k,Iter_MixingId,SubspaceInds,DistFunc,DegFunc,FitFunc,n,A);
% AllAhat=FnSubSpaceFind_MixingNew(ComplementOrthofSubSpaces,Th_RansacMixing,k,Iter_MixingId,SubspaceInds,DistFunc,DegFunc,FitFunc,n,A);
          
% Error estimation 
[AhatNew,M]=nearest2(Ahat,(A));
BAS=0;
for e=1:n
    dd=acos(AhatNew(:,e)'* A(:,e));
    BAS=dd+BAS;
end
Err=abs(BAS)
AhatNew



% L = n; % Maximum number of iterations for the internal steepest descent loop
% sigma_off = 0.001;
%
% % Important Algorithm's Parameters
% sigma_decrease_factor = 0.5;
% if sigma_off>0
%     sigma_min = sigma_off*4;
% else
%     sigma_min = 0.00001;
% end
%%
% Ahat=ComplementOrthofSubSpaces2;
% %UBI
% v=FnColNormalizer(randn(m,1));
% Sigma=0.2;
% v=FnMixingCalc_FMN3(SubSpaces,m,sigma,sigma_min, sigma_decrease_factor, L);
% v=FnMixingCalc_FMN4(SubSpaces,m);
% FnRound_h_v(SubSpaces,v,Sigma)
% [Ahat,MinEV]=FnMixingCalc_Ehsan(ComplementOrthofSubSpaces,k,c,n);
% [Ahat,MinEV]=FnMixingCalc_Fast(ComplementOrthofSubSpaces,k,c,n);

% Th3=1e-10;
% 
% Ahat=[];
% C_Ahat=[];
% AllAhat=[];
% N2=0;
% Th1=1e-13;
% combsubset=[];
% g=nchoosek(n-1,k-1);
% Cg=nchoosek(c,g);
% Mtotal=1e-7*Cg;
% % %%
% while size(combsubset,1)<Cg %N<NMax
%  rcmbn=FnTest_random_nchoosek(c,g,Mtotal);
%  combsubset=[combsubset; rcmbn'];
%  combsubset=unique(combsubset,'rows');
%  size(combsubset,1)
% end
% [Ahat,MinEV]=FnMixingCalc_Recursive(ComplementOrthofSubSpaces,n,combsubset);



% Ahat=v;
% [Ahat,MinEV]=FnMixingCalc_Ehsan_V2(ComplementOrthofSubSpaces,k,c,n);
% [Ahat]=FnCS2_New(ComplementOrthofSubSpaces,m,Th3);
% % [indd,Ahat]=kmeans(Ahat',n);
% % Ahat=Ahat';
% C_Ahat=[];
% Th1=1e-15;
% NN=1;
% AllAhat=[];
% for j=1:size(Ahat,2)
%     [C_Ahat,Winner]=FnBBC(Ahat(:,j),C_Ahat,NN,Th1);%Alg2:BBC
%     NN=NN+1;
%     ChannelNum=size(C_Ahat,2)
%
%     AllAhat=[AllAhat C_Ahat];
% end
% [Ahat,MinEV]=FnMixingCalc_BM(ComplementOrthofSubSpaces,k,c,n,QRSubspaceInds,A);
%  [Ahat,MinEV]=FnMixingCalc_BM_EE2(ComplementOrthofSubSpaces,k,c,n,QRSubspaceInds,A);


% figure, plot(sort(MinEV,'descend'),'*')
% [Ahat, MinEV]=FnMixingCalc_BM(ComplementOrthofSubSpaces,k,c,n,SubspaceInds,A);
% A=H;

