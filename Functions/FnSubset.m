function Out=FnSubset(n, k)
if k == 0
    Out= 1;
end
    if n == k
        Out= 1;
    else
       
        Out= nchoosek(n-1, k-1) + nchoosek(n-1, k);
    end