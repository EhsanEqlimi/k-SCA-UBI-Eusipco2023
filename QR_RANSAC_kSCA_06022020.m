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
m=3; % The number of sensors a.k.a the ambient dimension
n=4; % The number of sources
k=m-1; % The number of acitive source in each time point (nonzero elements
%in each col of S) or the dimension of each subspace
r=m-k; % The dimension of the normal space (normal vector/orthogonal complement) of each subspace
c=nchoosek(n,k); % The number of all possible r-dim subspaces
T= 1500; % The Number of data poins a.k.a the number of samples (time points)
Sigma=1e-9; % Parameter controls the standard deviation of normal noise over the
% zero sources (non-strict saprsity)
AMode=1; % k-SCA Condition for A are satisfied
IterNum=500;% The number of iteration to generate a good mixing matrix
RankTh=0.1; % to generate a good mixing matrix
Thr=1e-6;% RANSAC Threshold in Distance Function
Th3=1e-8;
Sigma1=0.02;
%% RANSAC Functions:
DistFunc = @FnDistanceClacBetweenSubspaceVector;
DegFunc= @(x) FnDegenerateEval(x,k);
g=nchoosek(n-1,k-1);
DegFunc2= @(x) FnDegenerateEval(x,g-1);

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
% figure, scatter(X(1,:),X(2,:));
% X=X+0.001*randn(m,T);
% [Sjj]=L1_norm_min(X,A);
%%
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
%% RANSAC based k-SCA
% Test1: k-SCA mixing mode (k=max i.e. m-1)
% load Sig.mat
% X=reshape(Frame_X_Fs,[512*65,2]);
% X=X';
% k=1;
% Subspace finding
[QRSubspaceInds,Clusters,SubSpaces,ComplementOrthofSubSpaces,ConnMat]=FnSubSpaceFind(X,Thr,k,c,SubspaceInds,DistFunc,DegFunc,FitFunc);

[QRSubspaceInds2,Clusters2,SubSpaces2,ComplementOrthofSubSpaces2,ConnMat2]=FnSubSpaceFind_Mixing(ComplementOrthofSubSpaces,1e-6,k,100,SubspaceInds,DistFunc,DegFunc,FitFunc,A);

% [Ahat,MinEV]=FnMixingCalc_Ehsan(ComplementOrthofSubSpaces,k,c,n);
Ahat=ComplementOrthofSubSpaces2;
%UBI
v=FnColNormalizer(randn(m,1));
Sigma=0.2;
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
[AhatNew,M]=nearest2(Ahat,(A));
BAS=0;
for e=1:n
    dd=acos(AhatNew(:,e)'* A(:,e));
    BAS=dd+BAS;
end
Err=abs(BAS)

AhatNew
Clusters(find(Clusters==0))=1;
[EstLabels] = bestMap(Labels,Clusters);
ClusterinError = (sum(Labels(:) ~= EstLabels(:)) / length(Labels))*100


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

