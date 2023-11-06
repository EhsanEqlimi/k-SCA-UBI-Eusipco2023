function [min_krank,mineig,min_sel,minvec,all_eig,all_sel,all_vec,worst_case,krank_old,mineig_old,min_sel_old,minvec_old]=k_rank_EE(A,thr,maxlevel)
n=size(A,2);
m=size(A,1);

if(~exist('maxlevel'))
   maxlevel=min(n,m); 
end

%min_krank=maxlevel;
min_krank=min(m,n);
worst_case=[];

if(~exist('thr'))
thr=1e-3;
end

mineig=1e20;
minvec=[];
min_sel=[];
krank=1;
st=1;
for i=1:min(n,min(m,maxlevel))%min(n,m):-1:2

           if(st)
               ii=min(n,min(m,maxlevel));
%          out=combntns(1:n,ii); 
           out=nchoosek(1:n,ii);
         j=1;
                svx=mnorm(svd(A(:,out(j,:)))) ;

       all_eig(:,j,ii)=mineig;    
       
       all_vec(:,j,ii)=svx; 
       all_sel(:,j,ii)=out(j,:);    
       st=0;
           end
    
        
%if(m<=n)
% out=combntns(1:n,i);
out=nchoosek(1:n,i);
%end
    for j=1:size(out,1)
       svx=mnorm(svd(A(:,out(j,:)))) ;
       if(svx(end)<mineig)
           
          mineig_old=mineig;
          minvec_old=minvec;
          min_sel_old=min_sel;
          krank_old=krank;
           
          mineig=svx(end);
          minvec=svx;
          min_sel=out(j,:);
          krank=i;
          
       end
       
       all_sel(1:size(out,2),j,i)=out(j,:);    
       all_vec(1:length(svx),j,i)=svx;
       all_eig(:,j,i)=svx(end);    
       
if(svx(end)<=thr),
    min_krank=min(i-1,min_krank);
    zx=zeros(1,maxlevel);
    zx(1:length(out(j,:)))=out(j,:);
    worst_case=[worst_case;zx];
    if(thr>0)
    break;
    end
    
end

           end
    
if(mineig<=thr),
   % worst_case=out;
   if(thr>0)
    break;
   end
end

end
%krank,mineig,min_sel,minvec
if 0&(~isempty(worst_case))
           mineig=mineig_old;
          minvec=minvec_old;
          min_sel=min_sel_old;
          krank=i-1;%krank_old;
end

if (~isempty(worst_case))
%krank=i-1;
end
    
    