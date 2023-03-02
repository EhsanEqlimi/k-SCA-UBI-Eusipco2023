function [Ahat,MinEV]=FnMixingCalc_BM_EE2(w,k,c,n,sub,A)
g=nchoosek(n-1,k-1);
% Cx=nchoosek(1:c,g);
% Cg=size(Cx,1);
% f=nchoosek(n,g);
% Cg=nchoosek(c,g);
AllPerms=nchoosek(1:size(w,2),size(w,1));
for i=1:size(AllPerms,1)%n%Cg
%     idx=find(prod((sub-i)')==0);
 wbar = w(:,AllPerms(i,:));
    %wbar=w(:,Cx(i,:));
%     wbar=w(:,Selected);
    R=wbar*wbar';
    [V,D]= eig(R);
    [Val,Ix] =sort(diag(D));
    MinEV(i)=abs(Val(1))/abs(Val(end));
    MinEV(i)=abs(Val(1));
    MinInd(i)=Ix(1);
    VV(:,i)=V(:,Ix(1));
end
[Val,Ix] =sort( MinEV);
f=n;
% temp=VV(:,Ix(1:Cg));
Ahat=VV(:,Ix(1:f));