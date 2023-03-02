function [C,Winner]=FnBBC(X,C,Itert,Th)
%%option1:
% X=fscolnorm(X);
%%%%%%%%%%%%%%%%%
if Itert ==1
    C(:,1)=X;
    Winner=1;
else
    %     for i=1:size(C,2)
    %         %option 2:
    %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %         [Distance(i),Ind(i)]=min([norm(C(:,i)-X) norm(C(:,i)+X)]);
    %     end
    %     [Value,Index]=min(Distance);
    %     if Ind(Index)==1
    %         X=X;
    %     else
    %         X=-X;
    %         %%%%%%%%%%%%%%%%%%%%%%%%%%
    %     end
    %     if Value <Th
    %         C(:,Index)=(C(:,Index)+X)./2;
    %         Winner= Index;
    %     else
    %         Winner=size(C,2)+1;
    %         C(:,size(C,2)+1)=X;
    %     end
    %% Option 3: Added by Ehsan 10/12/2014
 
    for i=1:size(C,2)
        CosineSimilarity=C(:,i)'*X;
        CosineSimilarity=CosineSimilarity/(norm(C(:,i))*norm(X));
        Dis(i)=1-abs(CosineSimilarity);
        Sgn(i)=sign(CosineSimilarity);
    end 
    [Value,Index]=min(Dis);
    if Value <Th
        C(:,Index)=(C(:,Index)+Sgn(Index)*X)./2;
        Winner= Index;
    else
        Winner=size(C,2)+1;
        C(:,size(C,2)+1)=Sgn(Index)*X;
    end
end

