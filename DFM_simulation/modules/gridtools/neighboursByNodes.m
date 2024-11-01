function N = neighboursByNodes(G)
%Derive neighbourship from common node (vertex) relationship
%
% SYNOPSIS:
%   N = neighboursByNodes(G)
%
% PARAMETERS:
%   G - Grid structure.
%
% RETURNS:
%   N - A neighbourship relation (m-by-2 array of inter-cell connections
%       (i.e., cell pairs)).  Two cells are defined to be neighbours (and
%       entered into the neighbourship relation) if they share a common
%       node.  The neighbourship relation is not symmetric and does not
%       include self-connections.
%
%       Use the statement
%
%          Adj = sparse([N(:,1); N(:,2); (1:G.cells.num).'], ...
%                       [N(:,2); N(:,1); (1:G.cells.num).'],  1 );
%
%       to convert the neighbourship relation 'N' into a (symmetric)
%       adjacency/Laplace matrix (neighbourhood), 'Adj'.
%
% NOTE:
%   This function uses SORTROWS and is potentially *very* expensive in
%   terms of memory use.  As an example, the statement
%
%      N = neighboursByNodes(cartGrid([60, 220, 85]))
%
%   requires about 6.5 GB of RAM and takes in the order of 20 seconds on a
%   workstation from the autumn of 2009.
%
% SEE ALSO:
%   sortrows.

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

   % Construction inspired by computeMultiPointTrans>createMapping.
   %
   c = rldecode(G.faces.neighbors, diff(G.faces.nodePos));
   n = repmat  (G.faces.nodes    , [2, 1]);

   i = c ~= 0;

   t = rlencode(sortrows(double([n(i), c(i)])));
   n = accumarray(t(:,1), 1);

   [i, j] = blockDiagIndex(n, n);
   N = sort([t(i,2), t(j,2)], 2);
   N = unique(N(N(:,1) ~= N(:,2), :), 'rows');
end
