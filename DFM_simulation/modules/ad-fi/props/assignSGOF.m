function f = assignSGOF(f, sgof, reg)
f.krG  = @(sg, varargin)krG(sg, sgof, reg, varargin{:});
f.krOG = @(so, varargin)krOG(so, sgof, f, reg, varargin{:});
f.pcOG = @(sg, varargin)pcOG(sg, sgof, reg, varargin{:});
end

function v = krG(sw, sgof, reg, varargin)
satinx = getRegMap(sw, reg.SATNUM, reg.SATINX, varargin{:});
T = cellfun(@(x)x(:,[1,2]), sgof, 'UniformOutput', false);
T = extendTab(T);
v = interpReg(T, sw, satinx);
end

function v = krOG(so, sgof, f, reg, varargin)
satinx = getRegMap(so, reg.SATNUM, reg.SATINX, varargin{:});
T = cellfun(@(x)x(:,[1,3]), sgof, 'UniformOutput', false);
T = extendTab(T);
if isfield(f, 'sWcon')
    v = interpReg(T, 1-so-f.sWcon, satinx);
else
    v = interpReg(T, 1-so, satinx);
end
end

function v = pcOG(sg, sgof, reg, varargin)
satinx = getRegMap(sg, reg.SATNUM, reg.SATINX, varargin{:});
T = cellfun(@(x)x(:,[1,4]), sgof, 'UniformOutput', false);
T = extendTab(T);
v = interpReg(T, sg, satinx);
end

