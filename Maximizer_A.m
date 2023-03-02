% This function performs part 2-c-optimization of the second part of the algorithm (mixing vector estimation) %%%
% Note that the optimization used in this function is the steepect ascent
% method (appendix 3 of the paper).
function [v,erhat]=Maximizer_A(v, miu, q, B, sigma);
erhat=0;
m=size(v,1);
N=size(B,3);
itr=0;
c=sigma^2;
while (itr<100 && c>=sigma^2/10)
    itr=itr+1;
    newV=v;
    grad=zeros(m,1);
    for(i=1:N)
        B1=B(:,:,i);
        D(:,i)=B1*(B1'*v)/sigma^2;
        landa(i)=exp(-(1-sum((B1'*v).^2))/(2*sigma^2));
    end;
    for(i=1:N)
        lan=landa;
        lan(i)=[];
        E=nsumk(lan,q-1);
        grad=grad+D(:,i)*landa(i)*E;
    end;
    if norm(grad)<1e-100
        erhat=1;
        return;
    end
    grad=grad/norm(grad);
    v=v+miu*grad;
    v=v/norm(v);
    c=min(norm(v-newV),norm(v+newV));
end