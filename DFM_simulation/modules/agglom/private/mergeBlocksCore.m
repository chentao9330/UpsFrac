function mrg = mergeBlocksCore(bN, bIVol, bIFlw, bINum, lbnd, ubnd, NB, cfac)
%Core implementation of MERGE primitive
%
% SYNOPSIS:
%   mrg = mergeBlocksCore(bN, bIVol, bIFlw, bINum, lbnd, ubnd, NB, cfac)
%
% PARAMETERS:
%   bN    - Coarse-scale neighbourship definition, possibly derived using
%           function 'blockNeighbourship'.
%
%   bIVol - Accumulated volume indicator per block.  If a particular block
%           is empty, its 'bIVol' value should be NaN.
%
%   bIFlw - Relative flow indicator per block.  If a particular block is
%           empty, its 'bIFlw' value should be NaN.
%
%   bINum - Number of cells in each block
%
%   lbnd  - Lower block volume bound.  This algorithm will merge a block B
%           into a neighbouring block if bIVol(B)<lbnd.
%
%   ubnd  - Upper block flow bound.  The merging process will attempt to
%           honour the criterion bIFlw(B)*bIVol(B)<=ubnd.
%
%   NB    - Upper bound on number of cells within a single block.
%
%   cfac  - When violating the soft constraints on block flow (and number
%           of cells within a block) by a factor cfac, the constraints(s)
%           become hard constraint(s).
%
% RETURNS:
%   mrg   - A merging operator.  Specifically, a NUMEL(bIVol)-element
%           vector containing (a subset of) the numbers 1:NUMEL(bIVol) such
%           that
%
%               p2 = mrg(p)
%
%           in which 'p' is a partition vector will create a new, possibly
%           uncompressed, partition vector 'p2' defining the merged blocks.
%           Use function 'compressPartition' to remove empty blocks from
%           this block defintion.
%
% SEE ALSO:
%   mergeBlocks2, blockNeighbourship, compressPartition.

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

   assert(cfac>=1, 'mergeBlocksCore:AssertionFailed',... 
      'Parameter <cfac> must be in the range [1,inf)');
   
   mrg = (1 : size(bIVol, 1)) .';

   conn = blockConnectivity(double(bN));

   for b = reshape(candidates(bIVol, lbnd), 1, []),
      if viable(b, conn, bIVol, bIFlw, lbnd),
         [mrg, conn, bIVol, bIFlw, bINum] = ...
            merge_into_neigh(b, mrg, conn, bIVol, bIFlw, bINum, ubnd, NB, cfac);
      end
   end
end

%--------------------------------------------------------------------------

function b = candidates(IVol, lbnd)
   b      = find(IVol < lbnd);
   [i, i] = sort(IVol(b));              %#ok
   b      = b(i);
end

%--------------------------------------------------------------------------

function tf = viable(b, conn, IVol, IFlw, lbnd)
   i  = blockNeighbours(conn, b);
   tf = (IVol(b) < lbnd) && any(isfinite(IFlw(i)));
end

%--------------------------------------------------------------------------

function [mrg, conn, IVol, IFlw, INum] = ...
      merge_into_neigh(b, mrg, conn, IVol, IFlw, INum, ubnd, NB, cfac)

   n   = blockNeighbours(conn, b);
   flw = IVol(b)*IFlw(b) + IVol(n).*IFlw(n);
   num = INum(b) + INum(n);

   feasible = ~( (flw > ubnd) | (num>NB) );

   if any(feasible),
      % Merge into feasible neighbour that most closely matches (block)
      % flow indicator of 'b'.

      n = n(feasible);
      [i, i] = min(abs(IFlw(b) - IFlw(n)));  %#ok
   else
      % Merge into neighbour that minimises violation of upper bounds.

      meas = flw/ubnd + num/NB;
      [i, i] = min(meas);                     %#ok
      if meas(i) > cfac, return, end;
   end

   into          = n(i);
   mrg(mrg == b) = into;

   % Update block indicator values for merging.
   % Recall:
   %   IVol is additive while IFlw is relative.
   IFlw(into) = IVol(into)*IFlw(into) + ...
                IVol( b  )*IFlw( b  );

   IVol(into) = IVol(into) + IVol(b)   ;   IVol(b) = inf;
   IFlw(into) = IFlw(into) / IVol(into);   IFlw(b) = inf;
   INum(into) = INum(into) + INum(b)   ;   INum(b) = inf;

   conn{into} = [conn{into}; conn{b}];
   conn{ b  } = [];
end
