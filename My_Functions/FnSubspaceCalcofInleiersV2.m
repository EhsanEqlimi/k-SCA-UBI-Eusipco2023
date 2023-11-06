function [Sub,OrthSub,MinSV]=FnSubspaceCalcofInleiersV2(X,Inliers,k)

[U,S,V]=svd(X(:,Inliers));
OrthSub=U(:,end);
Sub=U(:,1:end-1);
SS=S(find(S));
MinSV=min(SS);

