function [ClusterinError1,ClusterinError2]=FnSubspaceClusteringErrorFinder(Clusters,Labels)
Clusters(find(Clusters==0))=1;
[EstLabels] = bestMap(Labels,Clusters);
ClusterinError1 = (sum(Labels(:) ~= EstLabels(:)) / length(Labels))*100;

Ids=find(Clusters==0);
Clusters(Ids)=[];
Labels(Ids)=[];
[EstLabels] = bestMap(Labels,Clusters);
ClusterinError2 = (sum(Labels(:) ~= EstLabels(:)) / length(Labels))*100;