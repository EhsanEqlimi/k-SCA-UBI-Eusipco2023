%Farid Movahedi Naini UBI
function v=FnMixingCalc_FMN4(Sub,m)
% mu_0=100*Sigma;
Sigmas=[0.1,0.05,0.025,0.0125];
Mus=100*Sigmas.^2;
L=100;
for j=1:L
 v(:,j)=FnColNormalizer(randn(m,1)); 
 i=1;
g=FnGFunc(Sub,v(:,j),Sigmas(i));
 v(:,j)= v(:,j) + Mus(i)*g(j);
 v(:,j)=FnColNormalizer(v(:,j));
 


end