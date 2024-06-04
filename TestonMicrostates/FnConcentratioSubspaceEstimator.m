function B=FnConcentratioSubspaceEstimator(X,Sigma_B)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
% This part estimates the concentration subsapces (first part of the algorithm).
% Note that steps 1 and 2 are done in the file test.m. 
% Here steps 3 and 4 are combined together:
m=size(X,1);
n=size(X,2);
k=1;
N_B=58; %parameter N_B in the paper
L_B=10*N_B; %parameter L_B in the paper
L_A=15*n;   %parameter L_A in the paper   
q=4;    %parameter q in the paper  
std_noise=0.01;  %standard deviation of inactive sources
std_source=1;   %standard deviation of active sources
np=nchoosek(n,k);   %total number of concetration subspaces
T=30*np;    %length of the sources
TH1=.01;    %The treshold related to step 4 of the first part of the algorithm (subspace estimation)
TH2=0.03;   %The treshold related to step 3 of the second part of the algorithm (mixing vector estimation)
TH3=.1;     %The treshold related to step 4 of the second part of the algorithm (mixing vector estimation)
ndifB=0;
for j=1:L_B;
    j
    sub_B=randn(m,k);
    sub_B=orth(sub_B);
    for sigma=Sigma_B
        miu=10000*sigma^2;
        sub_B=Maximizer_B(X, sub_B, miu, sigma);
        sub_B=orth(sub_B);
    end

    if j==1
        B(:,:,1)=sub_B;
        ndifB=1;
    end

    flag=0;
    for i=1:ndifB
        R= sum(sqrt(sum((sub_B'-sub_B'*B(:,:,i)*B(:,:,i)').^2,2)))/k;
        if R<TH1
            flag=1;
        end   
    end
    if flag==0
        ndifB=ndifB+1;
        B(:,:,ndifB)=sub_B;
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This part does step 5 of the first part of the algorithm:

for i=1:length(B(1,1,:))
    costf(i)=CostfunctionF(X,B(:,:,i), sigma);
end
[y,index]=sort(-costf);
if length(B(1,1,:))>N_B
    B=B(:,:,index(1:N_B));
end