function PNV=FnPNVClac(X)

GSX=FnGramSchmidtOrth(X);
% [Q,R]=FnGramSchmidth_QR(X);
% Pesudo Normal vectors / Pesudo Complement Orthogonals !!!!!
I=eye(size(X,1));
PNV=I-GSX*GSX';