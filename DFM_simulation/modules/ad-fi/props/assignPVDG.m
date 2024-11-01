function f = assignPVDG(f, pvdg, reg)
f.bG  = @(pg, varargin)bG(pg, pvdg, reg, varargin{:});
f.BG  = @(pg, varargin)BG(pg, pvdg, reg, varargin{:});
f.muG = @(pg, varargin)muG(pg, pvdg, reg, varargin{:});
end

function v = bG(pg, pvdg, reg, varargin)
satinx = getRegMap(pg, reg.SATNUM, reg.SATINX, varargin{:});
T = cellfun(@(x)[x(:,1) 1./x(:,2)], pvdg, 'UniformOutput', false);
v = interpReg(T, pg, satinx);
end

function v = BG(pg, pvdg, reg, varargin)
satinx = getRegMap(pg, reg.SATNUM, reg.SATINX, varargin{:});
T = cellfun(@(x)x(:,[1,2]), pvdg, 'UniformOutput', false);
v = interpReg(T, pg, satinx);
end

function v = muG(pg, pvdg, reg, varargin)
satinx = getRegMap(pg, reg.SATNUM, reg.SATINX, varargin{:});
T = cellfun(@(x)x(:,[1,3]), pvdg, 'UniformOutput', false);
v = interpReg(T, pg, satinx);
end

