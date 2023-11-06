function [inliers,M,dist2]=FnDistanceClacBetweenSubspaceVector(P,x,t)
% if size(P,2)==1
%     P=P';
% end
if ~isempty(P)
    dist2=sum((P*x).^2);
    inliers=find(dist2<t);
    M=P;
else
    inliers=[];
    M=[];
end
