function [p,k]=FnSparseness(X)
for j=1:size(X,2)
    x=X(:,j);
    n=length(x);
    p(j)=1/(sqrt(n)-1)^-1*sqrt(n)-(sum(abs(x)))/sqrt(sum(x.^2));
    k(j)=p(j)*size(X,1);
end
