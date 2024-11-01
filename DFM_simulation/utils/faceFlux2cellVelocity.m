function v = faceFlux2cellVelocity(G, faceFlux)
%Transform face-based flux field to one constant velocity per cell.
%
% SYNOPSIS:
%   veclocity = faceFlux2cellVelocity(G, faceFlux);
%
% PARAMETERS:
%   G        - Grid structure.
%
%   faceFlux - Vector of fluxes corresponding to face ordering.
%
% RETURNS:
%   v        - d x G.cells.num matrix of cell velocities.
%
% SEE ALSO:
%   cellFlux2faceFlux, faceFlux2cellFlux.

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

%% Compute constant part of velocity field in each cell
   % Need cellNo
   cellNo = rldecode(1:G.cells.num, diff(G.cells.facePos), 2)';
   if ~isfield(G.cells, 'centroids'),
      G = computeGeometry(G);
   end
   C    = G.faces.centroids(G.cells.faces(:,1), :)-G.cells.centroids(cellNo, :);
   cf   = faceFlux2cellFlux(G, faceFlux);
   tmp  = bsxfun(@times, C, cf);
   v    = zeros(G.cells.num, size(G.nodes.coords, 2));
   for d = 1:size(G.nodes.coords, 2),
      v(:,d) = accumarray(cellNo, tmp(:,d))./G.cells.volumes;
   end
end
