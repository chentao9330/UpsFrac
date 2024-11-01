function p = compressPartition(p)
%Renumber coarse block partitioning to remove any empty coarse blocks.
%
% SYNOPSIS:
%   p = compressPartition(p)
%
% PARAMETERS:
%   p - Original partition vector, may contain empty coarse blocks.
%
% RETURNS:
%   p - Updated partition vector.
%       Renumbered so as to remove empty coarse blocks.
%
% NOTE:
%   If the original partition does not contain any empty coarse blocks,
%   applying this function is an expensive way of doing nothing.
%
% EXAMPLE:
%   p = partitionCartGrid([4, 4, 1], [2, 2, 2]);
%   [p, compressPartition(p)]
%
% SEE ALSO:
%   partitionCartGrid, partitionUI, processPartition, generateCoarseGrid.

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

   p = reshape(p, [], 1);

   if any(p < 1),
      error(msgid('Partition:NonPositive'), ...
           ['Partition vector ''p'' must contain positive integer', ...
            ' entries exclusively']);
   end

   active        = find(accumarray(p, 1) > 0);
   compr         = zeros([max(p), 1]);
   compr(active) = 1 : numel(active);

   p = compr(p);
end
