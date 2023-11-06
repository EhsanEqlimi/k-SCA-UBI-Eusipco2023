 function combsusbet=rand_nchoosek(N,g,samplesize,Mtotal,itrmax)

Nelements=N;
%Mtotal=12000
M=samplesize;
combsubset=[];
while(size(combsubset,1)<Mtotal) &&(itrmax>0)
    %combsubset=[combsubset;randi(Nelements,[M g])];
    %combsubset=[combsubset;sort(randi(Nelements,[M g])')'];
    dd=sort(randi(Nelements,[M g])');
    ddd=diff(dd);
    dd(:,find(prod(ddd)==0))=[];
    combsubset=[combsubset;dd'];
    combsubset=unique(combsubset,'rows');
    itrmax=itrmax-1;
end
combsusbet=combsubset(1:Mtotal,:);

