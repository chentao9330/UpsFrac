function bN = blockNeighbourship(N, p)
%Derive coarse-scale neighbourship from fine-scale information
%
% SYNOPSIS:
%   bN = blockNeighbourship(N, p)
%
% PARAMETERS:
%   N  - Fine-scale neighbourship definition.  Often, but not necessarily
%        always, equal to the 'faces.neighbors' field of a grid structure.
%
%   p  - Partition vector defining coarse blocks.
%
% RETURNS:
%   bN - Coarse-scale (block) neighbourship definition (unique inter-block
%        connections).
%
% NOTE:
%   This function uses SORTROWS.
%
% SEE ALSO:
%   grid_structure, generateCoarseGrid, sortrows.

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

   p  = [0; p];
   bN = p(N + 1);

   % Include inter-block connections whilst excluding boundary.
   bN = bN((bN(:,1) ~= bN(:,2)) & ~any(bN == 0, 2), :);

   % Extract unique inter-block connections (symmetry handled elsewhere).
   bN = unique(sort(bN, 2), 'rows');
end
