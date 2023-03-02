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
m=4; % The number of sensors a.k.a the ambient dimension
n=5; % The number of sources
k=3; % The number of acitive source in each time point (nonzero elements
m_Vec=[3,4,5,6,7];

%in each col of S) or the dimension of each subspace
r=m-k; % The dimension of the normal space (normal vector/orthogonal complement) of each subspace
c=nchoosek(n,k); % The number of all possible r-dim subspaces
T= c*30; %1000; % The Number of data poins a.k.a the number of samples (time points)
% Sigma=1e-3; % Parameter controls the standard deviation of normal noise over the
% zero sources (non-strict saprsity)
SigmaVec=[1e-5 1e-4 1e-3];
AMode=0; % k-SCA Condition for A is satisfied
IterNum=500;% The number of iteration to generate a good mixing matrix
RankTh=0.1; % to generate a good mixing matrix
% Th_RansacSubspace=1e-5 ;% RANSAC Threshold in Distance Function for subspace estimation
Th_RansacSubspaceVec=[1e-3 1e-4 1e-5 1e-6];  % RANSAC Threshold in Distance Function for subspace estimation

Th3=1e-8; % is not importanat!
ThBBC=1e-4; %BBC clustering for mixing idetification
Sigma1=0.02;% is not importanat!
Th_RansacMixing=1e-5;% RANSAC Threshold in Distance Function for mixing matrix identifiation
Iter_MixingId=50*c; % is not importanat!
p=0.95; % Desired probability that we get a good sample in ransac
e= 1-(1/c);%(c-1)/c;  % Probability that a point is an outlier
MaxIterEst=log(1-p)/log(1-(1-e)^k); %1e6;%log(1-p)/log(1-(1-e)^k);
if MaxIterEst<1000% Sometimes if k=1, it happnes
    MaxIterEst=1e4;
end
DualMode=0; %Twice subpspace clusteing
%% RANSAC Functions:
DistFunc = @FnDistanceClacBetweenSubspaceVector;%distance function
DegFunc= @(x) FnDegenerateEval(x,k);%For subsapce ransac
FitFunc=@(x) FnPNVClac(x);%RANSAC Model (Gram-schimdit)
g=nchoosek(n-1,k-1);

DegFunc2= @(x) FnDegenerateEval(x,m-1);% For mixing ransac
%% SCA Mixing (Simulted Sources)
%%'kSCA';%'kSCA';% Uniform Sparse Component Mixing  %'MSCA'=>Multiple Sparse Component Mixing, 'PermKSCA'
MixingMode='kSCANoisy';% {'kSCANoisy'%'PermkSCA'%kSCA';%'MSCA';%OnOffGauss};
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
%% Simple simulation
%generating the mixing matrix
% std_noise=0.01;  %standard deviation of inactive sources
% std_source=1;   %standard deviation of active sou
% A=randn(m,n);
% for j=1:n
%     A(:,j)=A(:,j)./norm(A(:,j));
% end
% %generating the sources
% p=k/n;
% z=rand(n,T);
% S=(z>p).*randn(n,T)*std_noise+(z<=p).*(randn(n,T)*std_source);
% for j=1:n
%     S(j,:)=S(j,:)-mean(S(j,:));
% end
% %generating the mixture matrix
% X=A*S;
%%
%% Sparse Component Mixing (Simulation)
for aa=1:length(SigmaVec)
[X,S,OrthA,A,Labels,SubspaceInds,SubspaceNum]=FnSparseComponentMixing(m,n,k,T,N,Nk,SigmaVec(aa),IterNum,RankTh,AMode,MixingMode);
%Omitting the points that are near the origin
g=sqrt(sum(X.^2,1));
X=X(:,g>.2);
Labels=Labels(g>.2);
%normalizing every column of the mixture matrix
R=sum(X.^2,1);
X=X./(ones(m,1)*sqrt(R));
% Ahat=PKDSC(X, n, k, Sigma_B, Sigma_A, N_B, q, L_B, L_A, TH1, TH2, TH3);
TrureConnMat=FnCreateConnMat(Labels);
figure('units','normalized','outerposition',[0 0 1 1])
% subplot(2,2,1)
% subplot_tight (2,2,1,[0.08 0.03])
% imagesc(TrureConnMat);title('True Connectivity Matrix');

% X=FnColNormalizer(X);% It is impotant when wh have clos orgin samples in X;
%% just test (it's not important)
% figure, scatter(X(1,:),X(2,:),X(2,:));
% X=X+0.001*randn(m,T);
% [Sjj]=L1_norm_min(X,A);
%% Subspace Estimation

for bb=1:length(Th_RansacSubspaceVec)
    
    
    [QRSubspaceInds,Clusters,SubSpaces,ComplementOrthofSubSpaces,ConnMat]=FnSubSpaceFind(X,Th_RansacSubspaceVec(bb),k,c,SubspaceInds,DistFunc,DegFunc,FitFunc,MaxIterEst,p);

if ~(length(find(Clusters==0))) && DualMode
    
    [QRSubspaceInds_step2,Clusters_step2,SubSpaces_step2,ComplementOrthofSubSpaces_step2,ConnMat_step2]=FnSubSpaceFind(X(:,find(Clusters==0)),Th_RansacSubspace,k,c,SubspaceInds,DistFunc,DegFunc,FitFunc,MaxIterEst,p);
    AllCOS=cat(numel(size(ComplementOrthofSubSpaces)),ComplementOrthofSubSpaces,ComplementOrthofSubSpaces_step2);
else
    AllCOS=   ComplementOrthofSubSpaces;
end
subplot(2,2,bb);
% subplot_tight (2,2,bb+1,[0.08 0.03])
[ClusterinError1,ClusterinError2]=FnSubspaceClusteringErrorFinder(Clusters,Labels);
imagesc(ConnMat);title(['Estimated Connectivity Matrix-' 'Th1=' num2str(Th_RansacSubspaceVec(bb)) sprintf('\n') 'Clusterin Error=' num2str(ClusterinError1) '%'])

end
end