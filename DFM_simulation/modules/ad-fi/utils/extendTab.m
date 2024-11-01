function T = extendTab(T, pm)
if nargin < 2
    pm = 1;
end
if iscell(T)
    for rn = 1:numel(T)
        t1 = T{rn}(1,:)  ; t1(1) = t1(1)-pm;
        te = T{rn}(end,:); te(1) = te(1)+pm;
        T{rn} = [t1; T{rn}; te];
    end
else
    t1 = T(1,:)  ; t1(1) = t1(1)-pm;
    te = T(end,:); te(1) = te(1)+pm;
    T = [t1; T; te];
end
            
