%Farid Movahedi Naini UBI
function g=FnMixingCalc_FMN(Sub,Sigma,n,m)
for j=1:100
 v=FnColNormalizer(randn(m,1));   
for i=1:size(Sub,3)
    d=FnDistanceMeasure(v',Sub(:,:,i));
    d2=d^2;
    Num=2*0.1^2;
    gtemp(i)=exp(-d2/Num);
end 
 g(j)=sum(gtemp) ;  
end
end