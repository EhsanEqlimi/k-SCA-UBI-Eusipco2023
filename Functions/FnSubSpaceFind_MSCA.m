function [QRSubspaceInds,Clusters,SubSpaces,ComplementOrthofSubSpaces,ConnMat]=FnSubSpaceFind_MSCA(X,Thr,k,c,SubspaceInds,DistFunc,DegFunc,FitFunc,MaxIterEst,p)
%"SubspaceInds" are not used to estimate the subspaces. Just for
%ecaluation!
SelectedIndices=1:size(X,2);
Clusters=zeros(1,size(X,2));
% Clusters=[];
m=size(X,1);
QRSubspaceInds=zeros(k,size(X,2));
ConnMat=zeros(size(X,2),size(X,2));
c1=nchoosek(n,k-1);
for i=1:c
   i
    % RANSAC
    [PNV, Inliers] = ransac(X(:,SelectedIndices), FitFunc, DistFunc, DegFunc, k, Thr,MaxIterEst,p);
    XSel=X(:,SelectedIndices);
    
    [Sub,OrthSub]=FnSubspaceCalcofInleiers(XSel,Inliers,k);
    %
%     Ans=XSel(:,Inliers)'*OrthSub;
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
    if ~iscell(SubspaceInds)
    QRSubspaceInds(:,SelectedIndices(Inliers))=repmat(SubspaceInds(i,:)',[1,length(SelectedIndices(Inliers))]);
    end
    SelectedIndices=setdiff(SelectedIndices,SelectedIndices(Inliers));
end

for i=1:c1
   i
    % RANSAC
    [PNV, Inliers] = ransac(X(:,SelectedIndices), FitFunc, DistFunc, DegFunc, k-1, Thr,MaxIterEst,p);
    XSel=X(:,SelectedIndices);
    
    [Sub,OrthSub]=FnSubspaceCalcofInleiers(XSel,Inliers,k);
    %
%     Ans=XSel(:,Inliers)'*OrthSub;
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
    if ~iscell(SubspaceInds)
    QRSubspaceInds(:,SelectedIndices(Inliers))=repmat(SubspaceInds(i,:)',[1,length(SelectedIndices(Inliers))]);
    end
    SelectedIndices=setdiff(SelectedIndices,SelectedIndices(Inliers));
end
