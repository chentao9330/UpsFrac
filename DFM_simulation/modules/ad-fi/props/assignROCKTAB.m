function f = assignROCKTAB(f, rocktab, reg)
f.pvMultR = @(p, varargin)pvMultR(p, rocktab, reg, varargin{:});
f.tranMultR = @(p, varargin)tranMultR(p, rocktab, reg, varargin{:});
end

function v = pvMultR(p, rocktab, reg, varargin)
rockinx = getRegMap(p, reg.ROCKNUM, reg.ROCKINX, varargin{:});
T = cellfun(@(x)x(:,[1,2]), rocktab, 'UniformOutput', false);
v = interpReg(T, p, rockinx);
end

function v = tranMultR(p, rocktab, reg, varargin)
rockinx = getRegMap(p, reg.ROCKNUM, reg.ROCKINX, varargin{:});
T = cellfun(@(x)x(:,[1,3]), rocktab, 'UniformOutput', false);
v = interpReg(T, p, rockinx);
end


