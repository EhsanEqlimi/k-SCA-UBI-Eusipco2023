% function [M,D1,D2,B1,B2] = distinctLRSCA(r)
% 
% Given r, this code generates a matrix M with
% (1) two distinct low-rank sparse decomposition of the form M = D_iB_i, 
%     where D_i in R^{r x r} and B_i in R^{r x (r^3-2r^2)} for i = 1,2, 
%     with ||B_i(:,j)||_0 = r-1 for all i,j, and  
% (2) r^2-2r points on each subspace spanned by r-1 columns of D_i (i=1,2)
%     that have spark r. 
%
% See Algorithm 1 in 
% J.E. Cohen and N. Gillis, 'Identifiability of Low-Rank Sparse Component 
% Analysis', August 2018. 

function [M,D1,D2,B1,B2] = distinctLRSCA_EE(m,n)

k = n-1; 
% c=nchoosek(n,k);
% Randomly picking normal vectors that will correspond to the hyperplanes
% spanned by r-1 columns of D1 and D2: 
Z1 = randn(m,n); 
Z2 = randn(m,n); 
disp('*****************'); 
% fprintf(' r = %1.0f \n', r); 
disp('*****************'); 
D1 = []; 
D2 = []; 
% Construct D1 and D2, as the intersection of r-1 hyperplanes
% whose normal vectors are given by columns Z1 and Z2, respectively. 
C = nchoosek([1:n],k); 
for i = 1 : size(C,1) 
    D1 = [D1 null( Z1(:,C(i,:))' )];
    D2 = [D2 null( Z2(:,C(i,:))' )];
end
% Construct points at the intersections of hyperplanes spanned by D1 and D2
M = [];  
B1 = []; 
B2 = []; 
r=n;
for i = 1 : r
    for j = 1 : r
        Z = [Z1(:,i) Z2(:,j)]; 
        % Find r-2 points in the intersection defined by { x | Z'*x = 0 }
        M_temp = null(Z'); 
        B1 = [B1 D1\M_temp]; 
        B2 = [B2 D2\M_temp];
        M = [M M_temp];
    end
end
% Check M = D1*B1 = D2*B2
if norm(D1*B1 - D2*B2,'fro') > 1e-9 || norm(M - D1*B1,'fro') > 1e-9 
    warning('Problem with the condition X = D1*B1 = D2*B2.') 
end
% Check sparisty of the columns of B1 and B2
if max( min( abs(B1) ) > 1e-9 ) || max( min( abs(B2) ) > 1e-9 )
    warning('B1 and B2 are not at least 1-sparse per column.'); 
end