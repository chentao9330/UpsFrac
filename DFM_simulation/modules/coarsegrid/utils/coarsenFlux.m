function cflux = coarsenFlux(cg, flux)
% Compute net flux on coarse faces
%
% SYNOPSIS:
%   cflux = coarsenFlux(cg, flux)
%
% PARAMETERS:
%   cg      - Coarse grid including parent grid.
%
%   flux    - Fine grid fluxes.
%
% RETURNS:
%   cflux   - Accumulated net flux on coarse faces
%
% REMARK:
%
% SEE ALSO:
%   coarsenGeometry, coarseConnections

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

   assert(isfield(cg, 'parent'), ...
      'Huh!? Field ''parent''missing in coarse grid.  Did you really supply a coarse grid?');
   
   assert(isnumeric(flux) && ...
          size(flux, 1) == cg.parent.faces.num && ...
          size(flux, 2)==1,...
          'Flux must be a numeric column vector of fine-grid face fluxes.');

   faceno = rldecode(1:cg.faces.num, diff(cg.faces.connPos), 2)';
   cflux  = accumarray(faceno, flux(cg.faces.fconn).*fineToCoarseSign(cg));
end
