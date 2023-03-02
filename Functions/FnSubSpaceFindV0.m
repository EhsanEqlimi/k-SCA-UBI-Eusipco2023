function [QRSubspaceInds,Clusters,SubSpaces,ComplementOrthofSubSpaces,ConnMat]=FnSubSpaceFindV0(X,Thr,k,c,SubspaceInds,DistFunc,DegFunc,FitFunc)
SelectedIndices=1:size(X,2);
Clusters=zeros(1,size(X,2));
% Clusters=[];
m=size(X,1);
QRSubspaceInds=zeros(k,size(X,2));
ConnMat=zeros(size(X,2),size(X,2));
for i=1:c
   i
    % RANSAC
    [PNV, Inliers] = ransacV0(X(:,SelectedIndices), FitFunc, DistFunc, DegFunc, k, Thr);
    XSel=X(:,SelectedIndices);
    
    [Sub,OrthSub]=FnSubspaceCalcofInleiers(XSel,Inliers,k);
    %
    Ans=XSel(:,Inliers)'*OrthSub;
    if k==1
        SubSpaces(:,i)=Sub;
        ComplementOrthofSubSpaces(:,:,i)=OrthSub;
    elseif m-k==1
        SubSpaces(:,:,i)=Sub;
        ComplementOrthofSubSpaces(:,i)=OrthSub;
    else
        SubSpaces(:,:,i)=Sub;
        ComplementOrthofSubSpaces(:,:,i)=OrthSub;
    end
   
    Clusters(SelectedIndices(Inliers))=i;
    ConnMat(SelectedIndices(Inliers),SelectedIndices(Inliers))=1;
    QRSubspaceInds(:,SelectedIndices(Inliers))=repmat(SubspaceInds(i,:)',[1,length(SelectedIndices(Inliers))]);
    SelectedIndices=setdiff(SelectedIndices,SelectedIndices(Inliers));
end

