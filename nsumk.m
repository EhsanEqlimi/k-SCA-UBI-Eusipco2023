% This function is a fast algorithm for computing function h (function
% number 5 in the paper).
% Note that the idea in the algorithm used in this function is the same as
% appendix 2 but its implementation is different to reduce the cost
% computation.
function su=nsumk(a,q);
for(r=1:q)
    a(r)=sum(a.^r);
end;
b(1)=a(1);
for(r=2:q-1)
    c(1,r)=a(1)*a(r)-a(r+1);
end;
if(q==1)
    su=a(1);
    return;
end;
b(2)=(a(1)^2-a(2))/2;
for(s=2:q-1)
    for(r=2:q-s)
        c(s,r)=b(s)*a(r)-c(s-1,r+1);        
    end;
    b(s+1)=(b(s)*a(1)-c(s-1,2))/(s+1);
end;
su=b(q);