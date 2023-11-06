function [Ahat,MinEV]=FnMixingCalc_EhsanNew(w,k,c,n,A)
g=nchoosek(n-1,k-1);
% Cx=nchoosek(1:c,g);%Slow
% Cx=VChooseK(1:c,g);%Fast (c code)
% Cg=size(Cx,1);
% f=nchoosek(n,g);
Cg=nchoosek(c,g);
samplesize=10; %findig combnts among samplesize number of random generated combnts
Mtotal=1e4;     % total number of unique found comntns at each trial
itrmax=1e5;    % to avoid infinty loop when we do not approach Mtotal


combsubset=[];
while(1)
    rcmbn=rand_nchoosek(c,g,samplesize,Mtotal,itrmax)';
    combsubset=[combsubset  rcmbn];
    combsubset=combsubset';
    for i=1:size(combsubset,1)
        if numel(size(w))>=3
            wbar=w(:,:,combsubset(i,:));
            wbar=reshape(wbar, [size(wbar,1),size(wbar,2)*size(wbar,3)]);
        else
            wbar=w(:,combsubset(i,:));
            %         wbar=w(:,Inds(1:g));
        end
        R=wbar*wbar';
        [V,D]= eig(R);
        [Val,Ix] =sort(diag(D));
        MinEV(i)=abs(Val(1))/abs(Val(end));
        MinEV(i)=abs(Val(1));
        MinInd(i)=Ix(1);
        VV(:,i)=V(:,Ix(1));
    end
    [Val,Ix2] =sort( MinEV);
    Sel=VV(:,Ix2(1));
end
end
% combsubset=unique(combsubset,'rows');
% size(combsubset,1)
% end
%
% for i=1:Cg
%     i
%     Inds=randperm(c);
%
%     if numel(size(w))>=3
%         %         wbar=w(:,:,Cx(i,:));
%         wbar=w(:,:,Inds(1:g));
%         wbar=reshape(wbar, [size(wbar,1),size(wbar,2)*size(wbar,3)]);
%     else
%         %         wbar=w(:,Cx(i,:));
%         wbar=w(:,Inds(1:g));
%     end
%     R=wbar*wbar';
%     [V,D]= eig(R);
%     [Val,Ix] =sort(diag(D));
%     MinEV(i)=abs(Val(1))/abs(Val(end));
%     MinEV(i)=abs(Val(1));
%     MinInd(i)=Ix(1);
%     VV(:,i)=V(:,Ix(1));
% end
% [Val,Ix] =sort( MinEV);
% f=n;
% % temp=VV(:,Ix(1:Cg));
% Ahat=VV(:,Ix(1:f));
