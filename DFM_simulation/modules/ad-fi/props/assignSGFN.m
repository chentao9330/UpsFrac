function f = assignSGFN(f, sgfn, reg)
f.krG  = @(sg, varargin)krG(sg, sgfn, reg, varargin{:});
f.pcOG = @(sg, varargin)pcOG(sg, sgfn, reg, varargin{:});
end

function v = krG(sg, sgfn, reg, varargin)
satinx = getRegMap(sg, reg.SATNUM, reg.SATINX, varargin{:});
T = cellfun(@(x)x(:,[1,2]), sgfn, 'UniformOutput', false);
T = extendTab(T);
v = interpReg(T, sg, satinx);
end

function v = pcOG(sg, sgfn, reg, varargin)
satinx = getRegMap(sg, reg.SATNUM, reg.SATINX, varargin{:});
T = cellfun(@(x)x(:,[1,3]), sgfn, 'UniformOutput', false);
T = extendTab(T);
v = interpReg(T, sg, satinx);
end

