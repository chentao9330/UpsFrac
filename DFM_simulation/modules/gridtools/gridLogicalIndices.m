function ijk = gridLogicalIndices(G, varargin)
% Given grid G and optional subset of cells, find logical indices.
%
% SYNOPSIS:
%   ijk = gridLogicalIndices(G)
%   ijk = gridLogicalIndices(G, c)
%
% PARAMETERS:
%   G    - Grid structure
%
% OPTIONAL PARAMETERS
%   c    - cells for which logical indices are to be computed. Defaults to
%          all cells (1:G.cells.num);
%
% RETURNS:
%   ijk  - cell array of size G.cartDims where ijk{1} contains logical
%          indices corresponding to G.cartDims(1) and so on.
%
% EXAMPLE:
%   G = cartGrid([2, 3]);
%   ijk = gridLogicalIndices(G, c);

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
    assert(isfield(G, 'cartDims'), ...
        ['Logical indices requires that the grid ', ...
        'has a valid ''cartDims'' field.']);

    if exist('narginchk', 'builtin'),
        narginchk(1, 2);
    else
        error(nargchk(1, 2, nargin, 'struct'));
    end
    
    if nargin == 2 && isnumeric(varargin{1}),
        cells = reshape(varargin{1}, [], 1);
    else
        cells = (1:G.cells.num) .';
    end
    ijk = cell(size(G.cartDims));

    [ijk{:}] = ind2sub(G.cartDims, G.cells.indexMap(cells));
end
