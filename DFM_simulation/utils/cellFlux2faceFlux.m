function faceFlux = cellFlux2faceFlux(G, cellFlux)
%Transform cell-based flux field to face-based.
%
% SYNOPSIS:
%   faceFlux = cellFlux2faceFlux(G, cellFlux)
%
% PARAMETERS:
%   G        - grid structure
%
%   cellFlux - Set of flux vectors in cell-wise ordering.  One column for
%              each flux vector.
%
% RETURNS:
%   faceFlux - Set of face-wise ordered fluxes.  One column for each of the
%              input flux vectors.
%
% SEE ALSO:
%   faceFlux2cellFlux.

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

cf       = G.cells.faces(:,1);
cellNo   = rldecode(1:G.cells.num, diff(G.cells.facePos), 2) .';
sgn      = 2*(G.faces.neighbors(cf, 1) == cellNo) - 1;

faceFlux = sparse(double(cf), 1:numel(cf), sgn) * cellFlux;
faceFlux = bsxfun(@rdivide, faceFlux, accumarray(cf, 1));
