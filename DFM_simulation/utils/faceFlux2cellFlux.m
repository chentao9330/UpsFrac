function cellFlux = faceFlux2cellFlux(G, faceFlux)
%Transform face-based flux field to cell-based.
%
% SYNOPSIS:
%   cellFlux = faceFlux2cellFlux(G, faceFlux);
%
% PARAMETERS:
%   G        - Grid structure.
%
%   faceFlux - Set of face-wise ordered fluxes.  One column for each flux
%              vector.
%
% RETURNS:
%   cellFlux - Set of flux vectors in cell-wise ordering.  One column for
%              each of the input flux vectors.
%
% SEE ALSO:
%   cellFlux2faceFlux.

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

cellNo   = rldecode(1:G.cells.num, diff(G.cells.facePos), 2)';
sgn      = 2*(G.faces.neighbors(G.cells.faces(:,1), 1) == cellNo) - 1;
cellFlux = bsxfun(@times, sgn, faceFlux(G.cells.faces(:,1), :));
