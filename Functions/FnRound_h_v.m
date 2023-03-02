function FnRound_h_v(Sub,v,Sigma)
Sigma2=Sigma^2;
for l=1:size(Sub,3)
   for i=1:size(Sub,2)
   b_i=Sub(:,i,l);
   g=FnGFunc(Sub,v,Sigma);
   Temp(:,i)=b_i*(b_i'*v)*g;
       
   end
   SUM(:,l)=sum(Temp,2);
end
SUM2=sum(SUM,2);
Term1=(1/Sigma2)*SUM2;
    