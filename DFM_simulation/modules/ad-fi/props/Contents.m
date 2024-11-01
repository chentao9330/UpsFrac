% Files
%  assignDENSITY    - dens of size ntpvtx3
%  assignPVCDO      - por  = pvcdo(pvtnum,1);  ref pres
%  assignPVDO       - f.muO = @(po, varargin)muO(po, pvdo, reg, varargin{:});
%  assignPVTW       - pwr  = pvtw(pvtnum,1);  ref pres
%  assignRelPerm    - if ~isfield(f, 'krOG')       two-phase water/oil
%  assignSOF3       - f.relperm3ph = @(sw, sg, varargin)relperm3ph(sw, sg, f, varargin);
%  initDeckADIFluid - props

%{
#COPYRIGHT#
%}
