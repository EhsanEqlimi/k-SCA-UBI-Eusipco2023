% This function calculates the cost function f (function 3 in the paper)
% for a given subspace and sigma and mixture matrix.
function res=CostfunctionF(X, B, sigma)
B=orth(B);
k=size(B,2);
y=-(ones(1,size(X,2))-sum((B'*X).^2,1))/(2*sigma^2);
res=sum(exp(y));