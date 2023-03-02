function [Sub,OrthSub]=FnSubspaceCalcofInleiers(X,Inliers,k)

[U,S,V]=svd(X(:,Inliers));
OrthSub=U(:,k+1:end);
Sub=U(:,1:k);

