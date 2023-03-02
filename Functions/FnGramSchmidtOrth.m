
function GSX=FnGramSchmidtOrth(X)
GSX=X(:,1)/norm(X(:,1));
I=eye(size(X,1));
for i=2:size(X,2)
    gsx=(I-GSX*GSX')*X(:,i);
    gsx=gsx/norm(gsx);
    GSX=[GSX gsx];
end
