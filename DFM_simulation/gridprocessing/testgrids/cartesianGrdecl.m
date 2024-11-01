function g = cartesianGrdecl(xi, yi, zi, varargin)
%Construct Cartesian grid with variable physical cell sizes.
%
% SYNOPSIS:
%   G = cartesianGrdecl(x, y)
%   G = cartesianGrdecl(x, y, 'depthz', dz)
%   G = cartesianGrdecl(x, y, z)
%   G = cartesianGrdecl(x, y, z, 'depthz', dz)
%
% PARAMETERS:
%   x,y,z    - Vectors giving cell vertices, in units of meters, of individual
%              coordinate directions.  Specifically, the grid cell at
%              logical location (I,J,K) will have a physical dimension of
%              [x(I+1)-x(I), y(J+1)-y(J), z(K+1)-z(K)] (meters).
%
%   dz       - Depth, in units of meters, at which upper reservoir nodes
%              are encountered.  Assumed to be a
%              NUMEL(x)-by-NUMEL(y) array of nodal depths.
%
%              OPTIONAL.
%              Default value: depthz = ZEROS([numel(x), numel(y)])
%                             (i.e., top of reservoir at zero depth).
%
%
% RETURNS:
%   grdecl - A GRDECL structure suitable for further processing by
%            function 'processGRDECL' or for permanent storage on disk by
%            function 'writeGRDECL'.
%
% EXAMPLE:
%   % 1) Create a 20-by-20-by-5 Cartesian grid and plot it.
%   grdecl = cartesianGrdecl(0:20, 0:20, 0:5);
%   plotGrid(processGRDECL(grdecl));
%   view(3), grid on, axis tight
%
%   % 2) Make a 20-by-20-by-5 Cartesian grid with wavy shape
%   [X, Y] = ndgrid(linspace(0, 2*pi/3, 21), linspace(0, 2*pi/5, 21));
%   grdecl = cartesianGrdecl(0:20, 0:20, 0:5, 'depthz', sin(X).*sin(Y));
%   plotGrid(processGRDECL(grdecl));
%   view(3), grid on, axis tight
%
% SEE ALSO:
%   processGRDECL, writeGRDECL, simpleGrdecl

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

% $Date: 2010-07-20 09:57:17 +0200 (Tue, 20 Jul 2010) $
% $Revision: 4989 $

   opt = struct('depthz', zeros(numel(xi), numel(yi)));
   opt = merge_options(opt, varargin{:});

   g.cartDims = [numel(xi), numel(yi), numel(zi)] - 1;
   [X, Y, Z]  = ndgrid(xi, yi, zi);
   Z = bsxfun(@plus, Z, reshape(opt.depthz, g.cartDims(1:2)+1));
   % Make pillars
   n     = prod(g.cartDims(1:2) + 1);
   lines = zeros([n, 6]);
   lines(:, [1, 4]) = reshape(X(:,:,[1, end]), [n, 2]);
   lines(:, [2, 5]) = reshape(Y(:,:,[1, end]), [n, 2]);
   lines(:, [3, 6]) = reshape(Z(:,:,[1, end]), [n, 2]);
   g.COORD = reshape(lines.', [], 1);

   % Assign z-coordinates
   % ind(d) == [1, 2, 2, 3, 3, ..., dims(d), dims(d), dims(d)+1]
   ind = @(d) 1 + fix((1 : 2*g.cartDims(d)) ./ 2);
   z   = Z(ind(1), ind(2), ind(3));


   g.ZCORN = z(:);
   g.ACTNUM = reshape(ones(g.cartDims, 'int32'), [], 1);
end
