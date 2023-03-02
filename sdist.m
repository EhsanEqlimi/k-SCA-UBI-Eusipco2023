% This function calculates the distance between two subspaces according to
% appendix 1.
function dis=sdist(B1, B2)
k=size(B1,2);
B1=orth(B1);
B2=orth(B2);
sum((B1'*B2).^2,2);
dis=sqrt(sum(ones(k,1)-sum((B1'*B2).^2,2)));