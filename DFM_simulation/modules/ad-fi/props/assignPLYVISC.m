function f = assignPLYVISC(f, plyvisc, reg)
f.muWMult = @(c, varargin)muWMult(c, plyvisc, reg, varargin{:});
end

function v = muWMult(c, plyvisc, reg, varargin)
satinx = getRegMap(c, reg.SATNUM, reg.SATINX, varargin{:});
plyvisc = extendTab(plyvisc);
v = interpReg(plyvisc, c, satinx);
end
