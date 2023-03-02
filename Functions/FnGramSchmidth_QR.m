function [Q,R]=FnGramSchmidth_QR(X)
Q=zeros(size(X));
R=zeros(size(X,2),size(X,2));
for i=1:size(X,2)
    gsx=X(:,i);
    for j=1:i-1
        R(j,i)=Q(:,j)'*X(:,i);
        gsx=gsx-R(j,i)*Q(:,j);
    end
    R(i,i)= norm(gsx);
    Q(:,i)=gsx/R(i,i);
            
end
    
    