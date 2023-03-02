function [QRSubspaceInds,Clusters,SubSpaces,ComplementOrthofSubSpaces,ConnMat]=FnSubSpaceFind_Mixing2(X,Thr,k,c,SubspaceInds,DistFunc,DegFunc,FitFunc,A)
% SelectedIndices=1:size(X,2);
SelectedIndices=randperm(size(X,2));
Clusters=zeros(1,size(X,2));
m=size(X,1);
QRSubspaceInds=zeros(k,size(X,2));
ConnMat=zeros(size(X,2),size(X,2));
C_Ahat=[];
N=0;
Th1=1e-7;
AllAhat=[];
E=0;
SavedOCSs=[];
It=0;
for i=1:c
   
   E=0;
    % RANSAC
    if length(SelectedIndices)>k
        E=E+1;
         ComplementOrthofSubSpaces=[];
        [PNV, Inliers] = ransac(X(:,SelectedIndices), FitFunc, DistFunc, DegFunc, k, Thr);
        XSel=X(:,SelectedIndices);
        
        [Sub,OrthSub]=FnSubspaceCalcofInleiers(XSel,Inliers,k);
        if k==1
            SubSpaces(:,i)=Sub;
            ComplementOrthofSubSpaces(:,:,E)=OrthSub;
        elseif m-k==1
            SubSpaces(:,:,i)=Sub;
            ComplementOrthofSubSpaces(:,E)=OrthSub;
        else
            SubSpaces(:,:,i)=Sub;
            ComplementOrthofSubSpaces(:,:,E)=OrthSub;
        end
        Clusters(SelectedIndices(Inliers))=E;
        ConnMat(SelectedIndices(Inliers),SelectedIndices(Inliers))=1;
        %     QRSubspaceInds(:,SelectedIndices(Inliers))=repmat(SubspaceInds(i,:)',[1,length(SelectedIndices(Inliers))]);
        SelectedIndices=setdiff(SelectedIndices,SelectedIndices(Inliers));
    else
        SavedOCSs=[SavedOCSs ComplementOrthofSubSpaces];
        SelectedIndices=randperm(size(X,2));
 
        for j=1:size(ComplementOrthofSubSpaces,2)
%             for t=1:size(ComplementOrthofSubSpaces,2)
            
            [C_Ahat,Winner(j)]=FnBBC2(ComplementOrthofSubSpaces(:,j),SavedOCSs,1e-13);%Alg2:BBC
%         end
%         
        
  end
       
        
%         AllAhat=[AllAhat C_Ahat];
        
 
    end
end

% N=2;
% ComplementOrthofSubSpaces( :, all( ~any( ComplementOrthofSubSpaces ), 1 ) ) =[];
% for j=1:size(ComplementOrthofSubSpaces,2)
%     if j==1
%         N=j
%     end
%     [C_Ahat,Winner(j)]=FnBBC(ComplementOrthofSubSpaces,C_Ahat,N,Th1);%Alg2:BBC
%     
%     ChannelNum=size(C_Ahat,2)
%     
%     AllAhat=[AllAhat C_Ahat];
%         end
QRSubspaceInds=[];