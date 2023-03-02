function g=FnGFunc(Sub,v,Sigma)

for f=1:size(Sub,3)
    d=FnDistanceMeasure(v',Sub(:,:,f));
    d2=d^2;
    Num=2*Sigma^2;
    gtemp(f)=exp(-d2/Num);
    
end 
 g=sum(gtemp) ; 