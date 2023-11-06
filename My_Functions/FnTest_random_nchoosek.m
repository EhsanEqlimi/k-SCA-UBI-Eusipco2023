function combsubset=FnTest_random_nchoosek(c,g,Mtotal)
% clear 
% N=10%35
% g=3%10
% Cx=nchoosek(N,g);

samplesize=10; %findig combnts among samplesize number of random generated combnts
% Mtotal=100;     % total number of unique found comntns at each trial
itrmax=1e4;    % to avoid infinty loop when we do not approach Mtotal


combsubset=[];
% while(1)
    rcmbn=rand_nchoosek(c,g,samplesize,Mtotal,itrmax)';
combsubset=[combsubset; rcmbn];
combsubset=unique(combsubset,'rows');
size(combsubset,1);
% end

