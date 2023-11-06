function PNV=FnPNVClac2(wbar)
R=wbar*wbar';
[V,D]= eig(R);
[Val,Ix] =sort(diag(D));
MinEV=abs(Val(1))/abs(Val(end));
MinEV=abs(Val(1));
MinInd=Ix(1);
PNV=V(:,Ix(1));
% GSX=FnGramSchmidtOrth(X);
% % [Q,R]=FnGramSchmidth_QR(X);
% % Pesudo Normal vectors / Pesudo Complement Orthogonals !!!!!
% I=eye(size(X,1));
% PNV=I-GSX*GSX';