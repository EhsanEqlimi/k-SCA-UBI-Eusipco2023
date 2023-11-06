function [Labels,SubspaceInd]=FnSubspaceLabelFinder(S,SubspaceInds)
for i=1:size(S,2)
   Ind=find(S(:,i));
   Ind=sort(Ind);
   Temp=abs(SubspaceInds-repmat(Ind',[size(SubspaceInds,1),1]));
   SubspaceInd{i}=SubspaceInds(find(sum(Temp,2)==0),:);
   Labels(i)=find(sum(Temp,2)==0);
   
end