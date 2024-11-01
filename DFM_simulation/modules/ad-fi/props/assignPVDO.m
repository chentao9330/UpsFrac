function f = assignPVDO(f, pvdo, reg)
f.BO  = @(po, varargin)BO(po, pvdo, reg, varargin{:});
f.bO  = @(po, varargin)bO(po, pvdo, reg, varargin{:});
%f.muO = @(po, varargin)muO(po, pvdo, reg, varargin{:});
f.BOxmuO = @(po, varargin) BOxmuO(po, pvdo, reg, varargin{:});
end

function v = BO(po, pvdo, reg, varargin)
pvtinx = getRegMap(po, reg.PVTNUM, reg.PVTINX, varargin{:});
T = cellfun(@(x)x(:,[1,2]), pvdo, 'UniformOutput', false);
v = interpReg(T, po, pvtinx);
end

function v = bO(po, pvdo, reg, varargin)
pvtinx = getRegMap(po, reg.PVTNUM, reg.PVTINX, varargin{:});
T = cellfun(@(x)[x(:,1), 1./x(:,2)], pvdo, 'UniformOutput', false);
v = interpReg(T, po, pvtinx);
end

function v = muO(po, pvdo, reg, varargin)
pvtinx = getRegMap(po, reg.PVTNUM, reg.PVTINX, varargin{:});
T = cellfun(@(x)x(:,[1,3]), pvdo, 'UniformOutput', false);
v = interpReg(T, po, pvtinx);
end

function v = BOxmuO(po, pvdo, reg, varargin)
  v1=BO(po, pvdo, reg, varargin{:});
  v2= muO(po, pvdo, reg, varargin);
  v=v1.*v2; 
end