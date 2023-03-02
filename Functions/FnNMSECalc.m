% MSE, Ehsan Eqlimi,01/07/2016,10/05/95
function [NMSE,NMSSum]=FnNMSECalc(A,Ahat)
N=min(size(A,2),size(Ahat,2));
for i=1:size(A,1)
    for j=1:N
        Num=sum((A(i,j)-Ahat(i,j)).^2);
%         if Num==0
%             Num=eps;
%         end
        Den=sum(A(i,j).^2);
        NMSE(i,j) = 10*log10(Num/Den);
    end
end
NMSSum=sum(NMSE(:));