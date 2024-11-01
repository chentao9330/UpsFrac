function plotContours(g, value, n)
% Plot contours of cell data.
%
% SYNOPSIS:
%   plotContours(g, value, n)
%
% DESCRIPTION:
%
% REQUIRED PARAMETERS:
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%
% RETURNS:
%
% NOTE:
%
% EXAMPLE
% 
% SEE ALSO:
%   plotCellData

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

levels = min(value): (max(value)-min(value))/n:max(value);
colors = round(1:size(colormap, 1)/n:size(colormap, 1));

for i=1:n,
   
   % fill space between contours
   I = [false; value >= levels(i) & value <= levels(i+1)];
   J = any(g.faces.neighbors==0, 2);
   f = find(any(I(g.faces.neighbors+1), 2) & J);
   plotFaces(g, f, colors(i), 'EdgeColor', 'none', 'facea', 0.4, 'outline', true);
   
   % Plot internal level sets
   I = [false; value <levels(i+1)];
   J = [false; value >levels(i+1)];
   f = find(any(I(g.faces.neighbors+1), 2) & any(J(g.faces.neighbors+1), 2));
   plotFaces(g, f, colors(i), 'EdgeColor', 'none', 'facea', 0.4, 'outline', true);
end

 
