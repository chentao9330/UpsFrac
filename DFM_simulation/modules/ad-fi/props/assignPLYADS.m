function f = assignPLYADS(f, plyads, reg)
f.ads = @(c, varargin)ads(c, plyads, reg, varargin{:});
end

function v = ads(c, plyads, reg, varargin)
satinx = getRegMap(c, reg.SATNUM, reg.SATINX, varargin{:});
plyads = extendTab(plyads);
v = interpReg(plyads, c, satinx);
end
