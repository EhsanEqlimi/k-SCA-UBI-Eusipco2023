function [QRSubspaceInds,Clusters,SubSpaces,ComplementOrthofSubSpaces,ConnMat,Ahat]=FnSubSpaceFind_Mixing5(X,ThBBC,Thr,k,c,SubspaceInds,DistFunc,DegFunc,FitFunc,n,A)
X=squeeze(X);
% if numel(size(X))>=3
%     X=reshape(X, [size(X,1),size(X,2)*size(X,3)]);
%
% end

Clusters=zeros(1,size(X,2));
m=size(X,1);
QRSubspaceInds=zeros(k,size(X,2));
ConnMat=zeros(size(X,2),size(X,2));
C_Ahat=[];
N=0;
% ThBBC=1e-3;
AllAhat=[];
E=0;
SavedCOS=[];
It=0;
% g=nchoosek(n-1,k-1);
p=0.95; % Desired probability that we get a good sample in ransac
e= 1-(1/c);%(c-1)/c;  % Probability that a point is an outlier
MaxIterEst= 1e7; %log(1-p)/log(1-(1-e)^4); %1e6;%log(1-p)/log(1-(1-e)^k);
SelectedIndices=1:size(X,numel(size(X)));
% for i=1:500;%c
i=0;
MaxN=0;
AllPNV=[];
while size(SavedCOS,2)<n && MaxN<3e3
   disp(['Mixing vector#'  num2str(size(SavedCOS,2)) '/' num2str(n)]) ;
% SavedCOS
%     E=E+1;
    % RANSAC
     MaxN=MaxN+1
    if length(SelectedIndices)>k  && size(SavedCOS,2)<n
       
        if numel(size(X))>=3
            Selected1=X(:,:,SelectedIndices);
        else
            Selected1=X(:,SelectedIndices);
        end
        
        Selected=reshape(Selected1, [size(Selected1,1),size(Selected1,2)*size(Selected1,3)]);
% %         Selected=reshape(Selected1, [size(Selected1,1)*size(Selected1,2),size(Selected1,3)]);

%         if numel(size(X))>=3
%             [PNV, Inliers] = ransacV0(Selected, FitFunc, DistFunc, DegFunc, m-1, Thr);
            [PNV, Inliers] = ransac(Selected, FitFunc, DistFunc, DegFunc, m-1 ,Thr,MaxIterEst,p);
            AllPNV=[AllPNV PNV];
           
%         else
%             [PNV, Inliers] = ransacV0(Selected, FitFunc, DistFunc, DegFunc, k, Thr);
%         end
%         if numel(size(X))>=3
%             XSel=X(:,:,SelectedIndices);
%         else
%             XSel=X(:,SelectedIndices);
%         end        %         [OrthSub,MinSV]=FnSubspaceCalcofInleiersV3(XSel,Inliers,k);
        if ~isempty(Inliers)
            %Inliers=Inliers(randperm(k));
                E=E+1;

            [Sub,OrthSub]=FnSubspaceCalcofInleiers(Selected,Inliers,m-1);

%             [Sub,OrthSub,MinSV(E)]=FnSubspaceCalcofInleiersV2(Selected,Inliers,m-1);
            %
            %             if k+1==1
            %                 SubSpaces(:,i)=Sub;
            %                 ComplementOrthofSubSpaces(:,:,E)=OrthSub;
            %             elseif m-k==1
            SubSpaces(:,:,E)=Sub;
            ComplementOrthofSubSpaces(:,E)=OrthSub;
            %             else
            %                 SubSpaces(:,:,i)=Sub;
            %                 ComplementOrthofSubSpaces(:,:,E)=OrthSub;
            %             end
            %             Clusters(SelectedIndices(Inliers))=E;
            %             ConnMat(SelectedIndices(Inliers),SelectedIndices(Inliers))=1;
            %     QRSubspaceInds(:,SelectedIndices(Inliers))=repmat(SubspaceInds(i,:)',[1,length(SelectedIndices(Inliers))]);
            %             SelectedIndices=setdiff(SelectedIndices,SelectedIndices(round(Inliers./size(X,2))));
%             SelectedIndices=randperm(size(X,numel(size(X))));
hh=randperm(15);
% RP=randperm(hh(1));
% SelectedIndices=randperm(size(X,numel(size(X)))-RP(1));
% if MaxN<700
%  SelectedIndices=randperm(size(X,numel(size(X))));
% elseif MaxN>700 && MaxN<1500 
% SelectedIndices=randperm(size(X,numel(size(X)))-5);
% elseif MaxN>1500 && MaxN<2000 
% SelectedIndices=randperm(size(X,numel(size(X)))-7);
% elseif MaxN>1500 && MaxN<2000 
% SelectedIndices=randperm(size(X,numel(size(X)))-4);
% elseif MaxN>2000
%    SelectedIndices=randperm(size(X,numel(size(X)))-3);
%  
% end

            for j=1:size(squeeze(ComplementOrthofSubSpaces),2)
                tempMAt=squeeze(ComplementOrthofSubSpaces);
                [SavedCOS,Winner]=FnBBC3(tempMAt(:,j),SavedCOS,ThBBC);%Alg2:BBC
                %            SavedCOS= [SavedCOS C_Ahat];
                
            end
            SavedCOS=squeeze(SavedCOS);
        end
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
            [SavedCOS,Winner]=FnBBC3(tempMAt(:,j),tempMAt,ThBBC);%Alg2:BBC
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