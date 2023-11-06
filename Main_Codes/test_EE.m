%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program is a sample file that tests the performance of the proposed
% algorithm of the extimation of the mixing matrix in multiple componnet SCA


clear all;
close all;
clc;

m=3;    %number of mixtures
n=5;   %number of sources
k=2;    %average number of active sources
% Sigma_B=[.15 .075 .037 0.018 ]; %The decreasing sequence of sigma in the first part of the algorithm (subspace estimation)
Sigma_B=[0.05 0.025 0.0125 ]; %The decreasing sequence of sigma in the first part of the algorithm (subspace estimation)

Sigma_A=[.1 .05 .025 .0125];    %The decreasing sequence of sigma in the second part of the algorithm (mixing vector estimation)
% N_B=58; %parameter N_B in the paper
N_B=10; %parameter N_B in the paper

% L_B=10*N_B; %parameter L_B in the paper
L_B=400; %parameter L_B in the paper
L_A=3;   %parameter L_A in the paper
q=4;    %parameter q in the paper
std_noise=0;  %standard deviation of inactive sources
std_source=1;   %standard deviation of active sources
np=nchoosek(n,k);   %total number of concetration subspaces
T= 2000;30*np;    %length of the sources
TH1=.01;    %The treshold related to step 4 of the first part of the algorithm (subspace estimation)
TH2=0.03;   %The treshold related to step 3 of the second part of the algorithm (mixing vector estimation)
TH3=.1;     %The treshold related to step 4 of the second part of the algorithm (mixing vector estimation)

%generating the mixing matrix
A=randn(m,n);
for j=1:n
    A(:,j)=A(:,j)./norm(A(:,j));
end

%generating the sources
p=k/n;
% p=k;
p=k/n;
z=rand(n,T);
S=(z>p).*randn(n,T)*std_noise+(z<=p).*(randn(n,T)*std_source);
% for j=1:n
%     S(j,:)=S(j,:)-mean(S(j,:));
% end

%generating the mixture matrix
X=A*S;

%omitting the points that are near the origin
g=sqrt(sum(X.^2,1));
X=X(:,g>.2);
%normalizing every column of the mixture matrix
R=sum(X.^2,1);
X=X./(ones(m,1)*sqrt(R));

% Estimating the mixing matrix
Ahat=PKDSC(X, n, k, Sigma_B, Sigma_A, N_B, q, L_B, L_A, TH1, TH2, TH3);
% error
[MixingIdentificatioError, MixingVectorerror,AhatNew] = FnMixingIdentificationError(A,Ahat) ;