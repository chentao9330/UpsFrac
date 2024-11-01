function f = assignPLYADS(f, plyvisc, reg)
f.muWMult = @(c, varargin)muWMult(c, plyvisc, reg, varargin{:});
end

function v = muWMult(c, plyvisc, reg, varargin)
pvtinx = getRegMap(c, reg.PVTNUM, reg.PVTINX, varargin{:});
v = interpReg(plyvisc, c, pvtinx);
end
