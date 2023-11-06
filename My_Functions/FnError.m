function Error=FnError(A,Ahat)
%% Error defined in Underdetermined mixing matrix estimation by exploiting sparsity of
%sources
for i=1:size(Ahat,2)
  Num(i)=norm(A(:,i)-Ahat(:,i)) ;
  Den(i)=norm(A)*norm(Ahat);
end

Error=sum(Num/Den);