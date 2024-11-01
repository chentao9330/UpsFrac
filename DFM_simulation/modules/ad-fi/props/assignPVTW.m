function f = assignPVTW(f, pvtw, reg)
ntpvt = numel(reg.PVTINX);
if ntpvt == 1
    f.cW  = pvtw(1, 3);
else
    f.cW  = pvtw(reg.PVTNUM, 3);
end
f.BW  = @(pw, varargin)BW(pw, pvtw, reg, varargin{:});
f.bW  = @(pw, varargin)bW(pw, pvtw, reg, varargin{:});
f.muW = @(pw, varargin)muW(pw, pvtw, reg, varargin{:});
end


function v = BW(pw, pvtw, reg, inx)
ntpvt = numel(reg.PVTINX);
if ntpvt == 1
    pvtnum = 1;
elseif nargin < 4
    pvtnum = reg.PVTNUM;
    assert(numel(pvtnum)==size(pw,1));
else
    pvtnum = reg.PVTNUM(inx);
end
pwr  = pvtw(pvtnum,1); % ref pres
bwr  = pvtw(pvtnum,2); % ref fvf
cw   = pvtw(pvtnum,3); % compress
X = cw.*(pw-pwr);
v = bwr./(1+X+X.^2/2);
end

function v = bW(pw, pvtw, reg, inx)
ntpvt = numel(reg.PVTINX);
if ntpvt == 1
    pvtnum = 1;
elseif nargin < 4
    pvtnum = reg.PVTNUM;
    assert(numel(pvtnum)==size(pw,1));
else
    pvtnum = reg.PVTNUM(inx);
end
pwr  = pvtw(pvtnum,1); % ref pres
bwr  = pvtw(pvtnum,2); % ref fvf
cw   = pvtw(pvtnum,3); % compress
X = cw.*(pw-pwr);
v = (1+X+X.^2/2)./bwr;
end

function v = muW(pw, pvtw, reg, inx)
ntpvt = numel(reg.PVTINX);
if ntpvt == 1
    pvtnum = 1;
elseif nargin < 4
    pvtnum = reg.PVTNUM;
    assert(numel(pvtnum)==size(po,1));
else
    pvtnum = reg.PVTNUM(inx);
end
pwr  = pvtw(pvtnum,1); % ref pres
muwr = pvtw(pvtnum,4); % ref visc
vbw  = pvtw(pvtnum,5); % viscosibility
Y = -vbw.*(pw-pwr);
v = muwr./(1+Y+Y.^2/2);
end

