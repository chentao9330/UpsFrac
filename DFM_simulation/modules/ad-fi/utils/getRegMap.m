function reginx = getRegMap(val, REGNUM, REGINX, cellInx)
nt = numel(REGINX);
if nt == 1
    reginx = {':'};
else
    if nargin < 4 %whole domain
        if size(val, 1)~=numel(REGNUM) %do not use numel in case of adi
            error('Region reference for input undefined');
        end
        reginx = REGINX; 
    else %reference to (small amount of) specific cells
        regnum = REGNUM(cellInx);
        if numel(cellInx) > 1
            if size(val, 1)~=numel(cellInx) %do not use numel in case of adi
                error('Number of cell indices must be same as input values');
            end
            reginx = cellfun(@(x)find(x==regnum), num2cell(1:nt), 'UniformOutput', false);
        elseif numel(cellInx) == 1 % allow single input (for exploring single cell functions)
            reginx = {ones(size(val))*regnum};
        else
            error('Got empty cellInx input. This is not happening...')    
        end
    end
end
end
