function h = plotFaceData(G, varargin)
%Plot face data on exterior grid faces to current axes (reversed Z axis).
%
% SYNOPSIS:
%       plotFaceData(G, data)
%       plotFaceData(G, data, 'pn1', pv1, ...)
%       plotFaceData(G, cells, data)
%       plotFaceData(G, cells, data, 'pn1', pv1, ...)
%   h = plotFaceData(...)
%
% PARAMETERS:
%   G       - Grid data structure.
%
%   data    - Vector of values for each face in G.
%
%   cells   - Vector of cell indices defining sub grid.  
%
%             If unspecified, function 'ploFaceData' will behave as if the
%             caller defined
%
%                 cells = 1 : G.cells.num
%
%             meaning graphical output will be produced for all cells in
%             the grid model 'G'.  If 'cells' is empty (i.e., if
%             ISEMPTY(cells)), then no graphical output will be produced.
%
%   'pn'/pv - List of property names/property values.  OPTIONAL.
%             This list will be passed directly on to function PATCH
%             meaning all properties supported by PATCH are valid.
%
% RETURNS:
%   h  - Handle to resulting patch object.  The patch object is added
%        directly to the current AXES object (GCA).
%        OPTIONAL.  Only returned if specifically requested.  If
%        ISEMPTY(cells), then h==-1.
%
% SEE ALSO:
%   plotCellData, plotFaces, patch, newplot.

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

   if mod(nargin, 2) == 0,
      cells    = 1:G.cells.num;
      facedata = varargin{1};
   else
      cells    = varargin{1};
      facedata = varargin{2};
   end
   
   if isempty(cells),
      warning(msgid('SubGrid:Empty'), ...
             ['Empty cell selection in ''plotFaceData''.', ...
              '  No graphics for you.'])
      if nargout > 0, h = -1; end
      return
   end
   
   assert(G.griddim == 3);

   f = boundaryFaces(G, cells);
   h = plotFaces(G, f, facedata(f), 'EdgeColor', 'none', ...
                 varargin{2 + mod(nargin, 2):end});
end
