function pf = cellPartitionToFacePartition(g, p, varargin)
% Construct partition of all faces in g from cell partition.
%
% SYNOPSIS
%    pf = cellPartitionToFacePartition(g, p, varargin)
%
% DESCRIPTION
%    Define unique partition id for each pair of cell partition ids
%    occurring in P, and construct a partitioning of all faces in g.
%
% PARAMETERS:
%   g       - Grid structure as described by grid_structure.
%
%   p       - Cell partition
%
% OPTIONAL PARAMETERS
%   'AllBoundaryFaces' - true/false determine if all fine-grid boundary
%                        faces should have a unique partition number to
%                        force them to be included (individually) in the
%                        coarse grid.  
%
% RETURNS:
%   pf      - Face partition, one integer per face in g.
%
% COMMENTS:
%
%
% SEE ALSO:
%   processFacePartition

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
   
   opt = struct('AllBoundaryFaces', false);
   opt = merge_options(opt, varargin{:});
   
   % Form face partition pf from cell partition p.
   p       = [0;p];
   [tmp,k] = sortrows(sort(p(g.faces.neighbors+1), 2));
   [n, n]  = rlencode(tmp);                                            %#ok
   pf      = zeros(g.faces.num, 1);
   pf(k)   = rldecode(1:numel(n), n, 2)';
   
   % Add all fine-faces on domain boundary
   if opt.AllBoundaryFaces,       
      ix     = any(g.faces.neighbors==0, 2);
      pf(ix) = max(pf)+(1:sum(ix))';
   end

end
