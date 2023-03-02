function [QRSubspaceInds,Clusters,SubSpaces,ComplementOrthofSubSpaces,ConnMat,Ahat]=FnSubSpaceFind_Mixing4(X,Thr,k,c,SubspaceInds,DistFunc,DegFunc,FitFunc,n,A)
% X=squeeze(X);
%     if numel(size(X))>=3
%         X=reshape(X, [size(X,1),size(X,2)*size(X,3)]);
% 
%     end
SelectedIndices=1:size(X,2);
Clusters=zeros(1,size(X,2));
m=size(X,1);
QRSubspaceInds=zeros(k,size(X,2));
ConnMat=zeros(size(X,2),size(X,2));
C_Ahat=[];
N=0;
Th1=1e-8;
AllAhat=[];
E=0;
SavedCOS=[];
It=0;
for i=1:c
    E=E+1;
    % RANSAC
    if length(SelectedIndices)>k
        [PNV, Inliers] = ransacV0(X(:,SelectedIndices), FitFunc, DistFunc, DegFunc, k, Thr);
        XSel=X(:,SelectedIndices);
%         [OrthSub,MinSV]=FnSubspaceCalcofInleiersV3(XSel,Inliers,k);
        [Sub,OrthSub,MinSV(E)]=FnSubspaceCalcofInleiersV2(XSel,Inliers,k);
        
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
        SavedCOS=squeeze(SavedCOS);
    elseif  size(SavedCOS,2)<n
        E=0;
        %         C_Ahat=ComplementOrthofSubSpaces;
        It=It+1;
        %         if numel(size(ComplementOrthofSubSpaces))>=3
        %
        %             ComplementOrthofSubSpaces_Temp=reshape(ComplementOrthofSubSpaces, [size(ComplementOrthofSubSpaces,1),size(ComplementOrthofSubSpaces,2)*size(ComplementOrthofSubSpaces,3)]);
        %         end
        for j=1:size(squeeze(ComplementOrthofSubSpaces),2)
            tempMAt=squeeze(ComplementOrthofSubSpaces);
            [SavedCOS,Winner]=FnBBC3(tempMAt(:,j),SavedCOS,Th1);%Alg2:BBC
            %            SavedCOS= [SavedCOS C_Ahat];
            
        end
        
        if isempty(SavedCOS)
            
            SavedCOS=  ComplementOrthofSubSpaces;
        end
        SelectedIndices=randperm(size(X,2));
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
Ahat=SavedCOS;