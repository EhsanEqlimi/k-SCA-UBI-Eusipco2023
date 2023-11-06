function [Sub,OrthSub,MinSV]=FnSubspaceCalcofInleiersV4(X,Inliers,k)

[U,S,V]=svd(X(:,Inliers));
OrthSub=U;
Sub=U(:,1:k);
SS=S(find(S));
MinSV=min(SS);

