function bc = psideh(bc, G, side, fluid, varargin)
% Impose hydrostatic Dirichlet boundary condition (pressure) on global side.
%
% SYNOPSIS:
%   bc = psideh(bc, G, side, fluid)
%   bc = psideh(bc, G, side, fluid, 'pn', pv)
%   bc = psideh(bc, G, side, fluid, I1, I2)
%   bc = psideh(bc, G, side, fluid, I1, I2, 'pn', pv)
%
% PARAMETERS:
%   bc     - boundary condition structure as defined by function 'addBC'.
%
%   G      - Grid structure as described by grid_structure.  Currently
%            restricted to grids produced by functions cartGrid and
%            tensorGrid.
%
%   side   - Global side from which to extract face indices.
%            Must be one of {'LEFT', 'RIGHT', 'FRONT', 'BACK', ...
%                            'BOTTOM', 'TOP'}
%
%   fluid      - fluid object to know densities
%
%   I1, I2 - Cell index ranges for local axes one and two, respectively.
%            Must contain *ALL* indices (e.g., I1 = 1:50, I2 = 1:4).
%
%            
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%   sat    - Fluid composition of fluid outside of the reservoir.
%            An m array of fluid phase saturations. If m=3, 'sat'
%            are interpreted as 1 <-> Aqua, 2 <-> Liquid, 3 <-> Vapor.
%
%            This field is to side the density of the outside fluid and to
%            set the saturation of incoming fluid in a transport solver
%
%            Default value: sat = 0 (assume single-phase flow).
%
%   range  - Restricts the search for outer faces to a subset of the cells
%            in the direction perpendicular to that of the face. Example:
%            if side='LEFT', one will only search for outer faces in the
%            cells with logical indexes [:,range,:].
%            Default value: range = [] (do not restrict search)
%
%   ref_depth - Reference depth for pressure. Default i 0
%
%   ref_pressure - Reference pressure. Default is 0
%
%
% RETURNS:
%   bc     - Updated boundary condition structure.
%
% EXAMPLE:
%   See simpleBC.m, simpleSRCandBC.m.
%
% SEE ALSO:
%   pside, fluxside, addBC, solveIncompFlow, grid_structure.

%{
Copyright 2009, 2010, 2011, 2012 SINTEF ICT, Applied Mathematics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

% $Date: 2012-12-19 16:47:24 +0100 (Wed, 19 Dec 2012) $
% $Revision: 10500 $

if ~isfield(G,'cartDims')
   error('psideh:NotImplemented', 'psideh is not implemented for this grid type');
end
cartDims = G.cartDims;
if length(cartDims)<3, cartDims(3)=1; end

error(nargchk(4, 14, nargin, 'struct'))

error(nargchk(4, 10, nargin, 'struct'))

if nargin==4 || ischar(varargin{1})
   switch upper(side),
      case {'TOP', 'BOTTOM'}
         I1 = 1:cartDims(1); I2 = 1:cartDims(2);
      case {'LEFT', 'RIGHT'}
         I1 = 1:cartDims(2); I2 = 1:cartDims(3);
      case {'FRONT','BACK'}
         I1 = 1:cartDims(1); I2 = 1:cartDims(3);
      otherwise,
         error(psideh:Side:Unknown',                  ...
            'Boundary side ''%s'' not supported in psideh', side);
   end
else
   I1 = varargin{1};
   I2 = varargin{2};
   varargin = varargin(3:end);
end

opt = struct('sat', 0, 'range', [],'ref_position',[0.0,0.0,0.0],'ref_pressure',0.0);
opt = merge_options(opt, varargin{:});
sat = opt.sat;

ix = boundaryFaceIndices(G, side, I1, I2, opt.range);

%assert (any(numel(pressure) == [1, numel(ix)]));
assert (numel(sat)==0 || any(size(sat,1) == [1, numel(ix)]));
assert( size(sat,1) == 1)
if size(sat,1)     == 1, sat      = sat(ones([numel(ix), 1]), :); end
%if numel(pressure) == 1, pressure = pressure(ones([numel(ix), 1])); end
%keyboard
%fluid.omega(struct('s',sat))+ ...
[mu,rho]=fluid.properties(struct('s',sat))
kr=fluid.relperm(sat);
lambda=bsxfun(@rdivide,kr,mu);
omega=sum(bsxfun(@times,lambda,rho),2)./sum(lambda,2);


pressure = sum( (G.faces.centroids(ix,:)- ...
                 repmat(opt.ref_position,length(ix),1)).*repmat(gravity(), ...
                 length(ix),1) ,2).*omega+...
                 opt.ref_pressure;

bc = addBC(bc, ix, 'pressure', pressure, 'sat', sat);
