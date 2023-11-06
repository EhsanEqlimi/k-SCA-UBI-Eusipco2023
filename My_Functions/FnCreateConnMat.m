function ConnMat=FnCreateConnMat(Labels)
ConnMat=zeros(max(Labels),max(Labels));
for i=unique(Labels)
    Ind=find(Labels==i);
    ConnMat(Ind,Ind)=1;
end