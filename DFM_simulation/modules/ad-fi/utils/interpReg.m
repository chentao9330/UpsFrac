function [yi, dyidxi] = interpReg(T, xi, reginx)
nreg = numel(reginx);
if nreg > 1
    yi = zeros(size(xi)); 
end
for k = 1:nreg
    if reginx{k} == ':' %for improved eff seperate this case
        yi = interpTable(T{k}(:,1), T{k}(:,2), xi);
    elseif ~isempty(reginx{k})
        yi(reginx{k}) = interpTable(T{k}(:,1), T{k}(:,2), xi(reginx{k}));
    end
end

if (nargout>1)
    dyidxi = dinterpReg(T, xi, reginx);
end
