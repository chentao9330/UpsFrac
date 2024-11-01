function resSol = initResSol(G, p0, s0, varargin)
%Initialize reservoir solution data structure.
%
% SYNOPSIS:
%   state = initResSol(G, p0, s0)
%
% DESCRIPTION:
%   Initialize the reservoir solution structure to uniform cell and face
%   pressures and no-flow flux fields throughout the reservoir.  Also, set
%   all-zero phase saturations and phase densities in all cells in model.
%
% PARAMETERS:
%   G  - Grid data structure.
%   p0 - Initial uniform reservoir pressure (scalar or a G.cells.num-by-1
%        vector).
%   s0 - Initial reservoir saturation (default: s0=0)
%
% RETURNS:
%   state - Initialized reservoir solution structure having fields
%              - pressure -- REPMAT(p0, [G.cells.num, 1])
%              - flux     -- Initial, all-zero face fluxes for all faces.
%              - s        -- REPMAT(s0, [G.cells.num, 1])
%
% REMARKS:
%   state.s   - An n-by-m array of fluid compositions with 'n' being the
%               number of faces in 'faces' and for m=3, the columns
%               interpreted as 1 <-> Aqua, 2 <-> Liquid, 3 <-> Vapor. This
%               field is for the benefit of transport solvers.
%
% SEE ALSO:
%   initWellSol, solveIncompFlow.

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

   if nargin <= 2, s0 = 0; end

   [nc, nf] = deal(G.cells.num, G.faces.num);

   if size(s0, 1) == 1,
      s0 = repmat(s0(1,:), [nc, 1]);
      
   elseif size(s0, 1) ~= nc,
      error('Initial saturation must either be [1 x np] or [G.clls.num x np].');
   end

   if numel(p0) == 1,
      p0 = repmat(p0, [nc, 1]);
   end
   
   resSol = struct('pressure', p0,             ...
                   'flux',     zeros([nf, 1]), ...
                   's',        s0);
   if nargin == 4, resSol.z = varargin{1}(ones([nc, 1]), :); end
end
