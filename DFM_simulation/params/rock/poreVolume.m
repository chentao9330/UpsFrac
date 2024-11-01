function pv = poreVolume(G, rock, varargin)
%Compute pore volumes of individual cells in grid.
%
% SYNOPSIS:
%   pv = poreVolume(G, rock)
%
% PARAMETER:
%   G    - Grid data structure.
%          Must contain valid field 'G.cells.volumes'.
%
%   rock - Rock data structure.
%          Must contain valid field 'rock.poro'.
%
% RETURNS:
%   pv   - Vector of size (G.cells.num)-by-1 of pore volumes for each
%          individual cell in the grid.  This typically amounts to the
%          expression rock.poro .* G.cells.volumes, but the function
%          handles non-unit net-to-gross factors as well.
%
% SEE ALSO:
%  computeGeometry.

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

pv = rock.poro .* G.cells.volumes;

if isfield(rock, 'ntg'),
   pv = pv .* rock.ntg;
end
