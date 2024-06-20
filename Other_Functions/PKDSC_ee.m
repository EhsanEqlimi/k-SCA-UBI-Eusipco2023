%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function estimates the mixing matrix. The argumants are defined in
% the file test.m and the related paper.

function [Ahat,B]=PKDSC_ee(X, n, k, Sigma_B, Sigma_A, N_B, q, L_B, L_A, TH1, TH2, TH3)
m=size(X,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
% This part estimates the concentration subsapces (first part of the algorithm).
% Note that steps 1 and 2 are done in the file test.m. 
% Here steps 3 and 4 are combined together:

ndifB=0;
for j=1:L_B;
%     j
Count=0;
    sub_B=randn(m,k);
    sub_B=orth(sub_B);
    for sigma=Sigma_B
        Count=Count+1;
        miu=10000*sigma^2;
        [sub_B,Dist_SubSub,f_B]=Maximizer_B_ee(X, sub_B, miu, sigma);
        
        sub_B=orth(sub_B);
        sub_B_All{j,Count}=sub_B;
        Dist_SubSub_All{j,Count}=Dist_SubSub;
        f_B_All{j,Count}=f_B;
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
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This part estimates the mixing vectors (second part of the algorithm).
% Note that step 1 is done in the file test.m. 
% Here steps 2 and 3 and 4 are combined together:

j=0;
Ahat=[];
while (size(Ahat,2)<n && j<L_A)
    j=j+1;
    %% start from a random vector
    i=0;
    erhat=0;
    v=randn(m,1);
    v=v/norm(v);
    while (erhat==0 && i<length(Sigma_A))           
        i=i+1;
        tempsigma=Sigma_A(i);
        step=tempsigma^2*100;
        [v,erhat]=Maximizer_A(v, step, q, B, tempsigma);
        v=v/norm(v);
%         erhat
    end
    if (erhat==0)
        for i=1:size(B,3)
            D_sub(i)=norm(v-B(:,:,i)*(B(:,:,i)'*v));
        end
        if sum(D_sub<TH2)>=q
            if size(Ahat,2)==0
                Ahat(:,1)=v;
            else
                for i=1:size(Ahat,2)
                    D_vec(i)=min(abs(acos(dot(v,Ahat(:,i)))),abs(acos(dot(v,-Ahat(:,i)))));
                end
                if min(D_vec)>TH3;
                    Ahat=[Ahat v];
                end
            end
        end
    end
end
                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       