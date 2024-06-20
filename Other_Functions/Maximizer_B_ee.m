% This function performs part 3-c-optimization of the first part of the algorithm (subspace estimation) %%%
% Note that the optimization used in this function is the steepect ascent
% method (appendix 3 of the paper).
function [B,Dist_SubSub,f_B]=Maximizer_B_ee(X, B, miu, sigma)
k=size(B,2);
m=size(B,1);
T=length(X(1,:));
itr=1;
c=sigma^2;
while(itr<101 && c>.01*sigma^2)
    old_B=B;
    fun=exp(-(1-sum((B'*X).^2))/(2*sigma^2));
    f_B(itr,:)=fun;
    for(i=1:k)
        grad(:,i)=-sum(((repmat(B(:,i)',[m 1])*X).*X).*(repmat(fun,[m 1])),2)/sigma^2;
    end;
    B=B+miu*grad/T;
    B=orth(B);
    itr=itr+1;
    c=abs(sdist(old_B,B));
    Dist_SubSub(itr)=c;
end