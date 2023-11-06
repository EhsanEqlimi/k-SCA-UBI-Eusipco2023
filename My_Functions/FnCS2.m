function [Ahat]=FnCS2(C,m,Ahat,Th)
F=0;
if numel(size(C))==3
    
    C=reshape(C,[size(C,1),size(C,2)*size(C,3)]);
end
if size(C,2)> m-1
    AllPerms=nchoosek(1:size(C,2),m);
    for i=1:size(AllPerms,1)
        Selected = C(:,AllPerms(i,:));
        [V,D]=eig(Selected*Selected');
        v=nchoosek(1:m,m-1);
        vv=(V(:,1));
        d=diag(D);
        d=d(2)/d(end);
        d3=max(max(abs(vv'*Selected)));
        for j=1:size(v,1)
            [V1,temp]=eig(Selected(:,v(j,:))*Selected(:,v(j,:))');
            temp=abs(diag(temp));
            d(j+1)=temp(2)/temp(end);
            vv(:,j+1)=V1(:,1);
            d3(j+1)=max(max(abs(V1(:,1)'*Selected)));
        end
        D=diag(D);
        D=abs(D);
        if abs(D(1)/D(end))<Th  && abs(D(2)/(D(end)+eps))>Th*1e4 %abs(D(2)/(D(1)+D(2)))>.999
            if max(d3)<Th %min(abs(d))>Th*1e8 && dd(end-1)/dd(end)<Th*1e4
                F=F+1;
                Ahat(:,F)=V(:,1);
                
            end
        end
    end
end