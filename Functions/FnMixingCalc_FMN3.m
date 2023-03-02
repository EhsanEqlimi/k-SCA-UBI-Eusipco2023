%Farid Movahedi Naini UBI
function v=FnMixingCalc_FMN3(Sub,m,sigma,sigma_min, sigma_decrease_factor, L)
% mu_0=100*Sigma;
for j=1:L
while sigma>sigma_min
    mu_0=100*sigma.^2;

 v(:,j)=FnColNormalizer(randn(m,1)); 
 N=1;
for i=1:size(Sub,3)
    d=FnDistanceMeasure(v(:,j)',Sub(:,:,i));
    d2=d^2;
    Num=2*0.1^2;
    gtemp(i)=exp(-d2/Num);
    
end 
 g(j)=sum(gtemp) ; 
 v(:,j)= v(:,j) + mu_0*g(j);
 v(:,j)=FnColNormalizer(v(:,j));
 
end
sigma = sigma * sigma_decrease_factor;

end