function [ComplementOrthofSubSpaces,Ahat]=FnSubSpaceFind_Mixing_Final(X,ThBBC,f,n,Thr,p,MaxIterEst)
X=squeeze(X);
m=size(X,1);% The number of sensors (ambient dimention).
c=size(X,numel(X)); % The number of subspaces (c=ncoosek(n,k)).
E=0;MaxN=0;
SavedCOS=[];
%AllPNV=[];
% p=0.95; % Desired probability that we get a good sample in ransac
% e= 1-(1/c);%(c-1)/c;  % Probability that a point is an outlier
% MaxIterEst= 1e7; %log(1-p)/log(1-(1-e)^4); %1e6;%log(1-p)/log(1-(1-e)^k);
% SelectedIndices=randperm(size(X,numel(size(X))));
SelectedIndices=[1:size(X,numel(size(X)))];

Level=m-1;
DistFunc = @FnDistanceClacBetweenSubspaceVector;%distance function
DegFunc= @(x) FnDegenerateEval(x,Level);%For subsapce ransac
FitFunc=@(x) FnPNVClac(x);%RANSAC Model (Gram-schimdit)
while size(SavedCOS,2)<n && MaxN<1500
%     disp(['Mixing vector#'  num2str(size(SavedCOS,2)) '/' num2str(n)]) ;
      MaxN=MaxN+1;
%       disp(MaxN);
    if numel(size(X))>=3
        Selected1=X(:,:,SelectedIndices);
    else
        Selected1=X(:,SelectedIndices);
    end
    Selected=reshape(Selected1, [size(Selected1,1),size(Selected1,2)*size(Selected1,3)]);
    %     PermutedInds=randperm(size(Selected,2));
    %     IndsF=PermutedInds(1:f);
    [~, Inliers] = ransac(Selected,FitFunc, DistFunc, DegFunc, Level ,Thr,MaxIterEst,p);
    if ~isempty(Inliers)
        
        [Sub,OrthSub]=FnSubspaceCalcofInleiers2(Selected,Inliers,m-1,1e-4);
        if ~isempty(Sub)
            E=E+1;
%             SubSpaces(:,:,E)=Sub;
            ComplementOrthofSubSpaces(:,E)=OrthSub;
            %         SelectedIndices=randperm(size(X,numel(size(X))));
            SelectedIndices=[1:size(X,numel(size(X)))];
            
            for j=1:size(squeeze(ComplementOrthofSubSpaces),2)
                tempMAt=squeeze(ComplementOrthofSubSpaces);
                [SavedCOS,~]=FnBBC3(tempMAt(:,j),SavedCOS,ThBBC);%Alg2:BBC in paper MSCA 2015, Eqlimi et al
            end
            SavedCOS=squeeze(SavedCOS);
        end
    end
end
Ahat=SavedCOS;