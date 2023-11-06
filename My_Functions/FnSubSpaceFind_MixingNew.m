function [AllAhat]=FnSubSpaceFind_Mixing(X,Thr,k,c,SubspaceInds,DistFunc,DegFunc,FitFunc,n,A)
SubSpaces=[];
g=nchoosek(n-1,k-1);
if numel(size(X))>=3
    
    X=reshape(X, [size(X,1),size(X,2)*size(X,3)]);
    
end
SelectedIndices=1:size(X,2);
Clusters=zeros(1,size(X,2));
m=size(X,1);
QRSubspaceInds=zeros(k,size(X,2));
ConnMat=zeros(size(X,2),size(X,2));
C_Ahat=[];
N=0;
Th1=1e-4;
AllAhat=[];
E=0;
SavedCOS=[];
It=0;
while size(AllAhat,2)<n
    [PNV, Inliers] = ransac(X(:,SelectedIndices), FitFunc, DistFunc, DegFunc, k, Thr);
    XSel=X(:,SelectedIndices);
    
    [Sub,OrthSub,MinSV]=FnSubspaceCalcofInleiers(XSel,Inliers,k);
%     ComplementOrthofSubSpaces=OrthSub;
%     
%     [AllAhat,Winner]=FnBBC4(ComplementOrthofSubSpaces,AllAhat,Th1);%Alg2:BBC
% %     AllAhat=[AllAhat SavedCOS];
    SelectedIndices=randperm(size(X,2));
    
end
%         %         if k==1
%         %             SubSpaces(:,i)=Sub;
%         %             ComplementOrthofSubSpaces(:,:,E)=OrthSub;
%         %         elseif m-k==1
%         %             SubSpaces(:,:,i)=Sub;
%         %             ComplementOrthofSubSpaces(:,E)=OrthSub;
%         %         else
%         %             SubSpaces(:,:,i)=Sub;
%         %             ComplementOrthofSubSpaces(:,:,E)=OrthSub;
%         %         end
%         Clusters(SelectedIndices(Inliers))=E;
%         ConnMat(SelectedIndices(Inliers),SelectedIndices(Inliers))=1;
%         %     QRSubspaceInds(:,SelectedIndices(Inliers))=repmat(SubspaceInds(i,:)',[1,length(SelectedIndices(Inliers))]);
%         SelectedIndices=setdiff(SelectedIndices,SelectedIndices(Inliers));
%     elseif  size(SavedCOS,2)<n
%         [Val,IDMin]=sort(MinSV);
%         ComplementOrthofSubSpaces=ComplementOrthofSubSpaces(:,IDMin(1));
%
%         ComplementOrthofSubSpaces=reshape(ComplementOrthofSubSpaces, [size(ComplementOrthofSubSpaces,1),size(ComplementOrthofSubSpaces,2)*size(ComplementOrthofSubSpaces,3)]);
%         ComplementOrthofSubSpaces=ComplementOrthofSubSpaces(:,end);
%         E=0;
%         %         C_Ahat=ComplementOrthofSubSpaces;
%         It=It+1;
%         %         if numel(size(ComplementOrthofSubSpaces))>=3
%         %
%         %             ComplementOrthofSubSpaces_Temp=reshape(ComplementOrthofSubSpaces, [size(ComplementOrthofSubSpaces,1),size(ComplementOrthofSubSpaces,2)*size(ComplementOrthofSubSpaces,3)]);
%         %         end
%         for j=1:size(ComplementOrthofSubSpaces,2)
%
%             [SavedCOS,Winner]=FnBBC3(ComplementOrthofSubSpaces(:,j),SavedCOS,Th1);%Alg2:BBC
%             %            SavedCOS= [SavedCOS C_Ahat];
%
%         end
%         if isempty(SavedCOS)
%
%             SavedCOS=  ComplementOrthofSubSpaces;
%         end
%         SelectedIndices=randperm(size(X,2));
%         %         AllAhat=[AllAhat C_Ahat];
%
%         ComplementOrthofSubSpaces=[];
%     end
% end
%
% % N=2;
% % ComplementOrthofSubSpaces( :, all( ~any( ComplementOrthofSubSpaces ), 1 ) ) =[];
% % for j=1:size(ComplementOrthofSubSpaces,2)
% %     if j==1
% %         N=j
% %     end
% %     [C_Ahat,Winner(j)]=FnBBC(ComplementOrthofSubSpaces,C_Ahat,N,Th1);%Alg2:BBC
% %
% %     ChannelNum=size(C_Ahat,2)
% %
% %     AllAhat=[AllAhat C_Ahat];
% %         end
% QRSubspaceInds=[];
% Ahat=SavedCOS;