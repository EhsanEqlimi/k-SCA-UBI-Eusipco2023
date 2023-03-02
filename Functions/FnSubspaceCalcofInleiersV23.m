function [OrthSub,MinSV]=FnSubspaceCalcofInleiersV23(X,Inliers,k)
%
% [U,S,V]=svd(X(:,Inliers));
% OrthSub=U(:,end);
% Sub=U(:,1:k);
% SS=S(find(S));
% MinSV=min(SS);
X=X(:,Inliers);
R=X*X';
[V,D]= eig(R);
[Val,Ix] =sort(diag(D));
MinEV=abs(Val(1))/abs(Val(end));
MinEV=abs(Val(1));
MinInd=Ix(1);
OrthSub=V(:,Ix(1));

