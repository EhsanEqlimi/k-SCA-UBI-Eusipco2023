%**************************************************************************
% Main file to run QR and RANSAC based k-SCA
% 20/12/94
% Ehsan Eqlimi,TUMS
%**************************************************************************
clc;
clear ;
close all;
warning off
%% Preliminaries
m=2; % The number of sensors a.k.a the ambient dimension
n=5; % The number of sources
k=m-1; % The number of acitive source in each time point (nonzero elements
%in each col of S) or the dimension of each subspace
r=m-k; % The dimension of the normal space (normal vector/orthogonal complement) of each subspace
c=nchoosek(n,k); % The number of all possible r-dim subspaces
T= 500; % The Number of data poins a.k.a the number of samples (time points)
Sigma=1e-5; % Parameter controls the standard deviation of normal noise over the
% zero sources (non-strict saprsity)
AMode=1; % k-SCA Condition for A are satisfied
IterNum=500;% The number of iteration to generate a good mixing matrix
RankTh=0.1; % to generate a good mixing matrix
Thr=1e-4;% RANSAC Threshold in Distance Function
%% RANSAC Functions:
DistFunc = @FnDistanceClacBetweenSubspaceVector;
DegFunc= @(x) FnDegenerateEval(x,k);
FitFunc=@(x) FnPNVClac(x);
%% SCA Mixing (Simulted Sources)
%%'kSCA';%'kSCA';% Uniform Sparse Component Mixing  %'MSCA'=>Multiple Sparse Component Mixing, 'PermKSCA'
MixingMode='PermkSCA'% {'kSCANoisy'%'PermkSCA'%kSCA';%'MSCA'};
Orth=0; %if Orth=1 A is orth
N=zeros(1,c);
N(1,:)=ceil(T/c); % N is the number of the subspaces in k-SCA mode
Nk=zeros(c,k);
for j=1:k
    for i=1:nchoosek(n,j)
        Nk(i,j)=ceil(ceil(T/k)/nchoosek(n,j));% N showes the number of each subspace in MSCA mode
    end
end

%% Sparse Component Mixing
[X,S,OrthA,A,Labels,SubspaceInds,SubspaceNum]=FnSparseComponentMixing(m,n,k,T,N,Nk,Sigma,IterNum,RankTh,AMode,MixingMode);
% X=X+0.001*randn(m,T);
% [Sjj]=L1_norm_min(X,A);
%% RANSAC based k-SCA
% Test1: k-SCA mixing mode (k=max i.e. m-1)
load Sig.mat
X=reshape(Frame_X_Fs,[512*65,2]);
X=X';
k=1;


[QRSubspaceInds,Clusters,SubSpaces,ComplementOrthofSubSpaces,ConnMat]=FnSubSpaceFind(X(:,1:10240),Thr,k,c,SubspaceInds,DistFunc,DegFunc,FitFunc);
% figure,imagesc(ConnMat);
Ahat=FnMixingCalc_Ehsan(ComplementOrthofSubSpaces,k,c,n);
% [Ahat, MinEV]=FnMixingCalc_BM(ComplementOrthofSubSpaces,k,c,n,SubspaceInds,A);
A=H;
[AhatNew,M]=nearest2(real(Ahat),A);
BAS=0;
for e=1:n
    dd=acos(AhatNew(:,e)'* A(:,e));
    BAS=dd+BAS;
end
Err=abs(BAS)
H
AhatNew
Clusters(find(Clusters==0))=1;
[EstLabels] = bestMap(Labels,Clusters);
ClusterinError = (sum(Labels(:) ~= EstLabels(:)) / length(Labels))*100
