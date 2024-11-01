function varargout = plotFaults(G, faults, varargin)
%Plot faults in model
%
% SYNOPSIS:
%       plotFaults(G, faults)
%   h = plotFaults(G, faults)
%
% PARAMETERS:
%   G      - Valid grid structure.
%
%   faults - Valid fault structure as defined by function 'processFaults'.
%
% RETURNS:
%   h - Two-element cell array of handles suitable for passing to GET or
%       SET.  Specifically, h{1} is a PATCH handle to the set of graphics
%       containing the fault faces while h{2} is a set of TEXT handles to
%       the corresponding fault names.
%
%       OPTIONAL.  Only returned if specifically requested.
%
% SEE ALSO:
%   plotFaces, patch, text.

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

   opt = struct('Height', 5, 'FontSize', 16, 'Color', 'r');
   opt = merge_options(opt, varargin{:});

   hf = plotFaces(G, vertcat(faults.faces), 'FaceColor', 'k');

   x  = cellfun(@(f) mean(G.faces.centroids(f,:), 1), ...
                { faults.faces }, 'UniformOutput', false);
   x  = vertcat(x{:});


   ht = text(x(:,1), x(:,2), x(:,3) - opt.Height,          ...
             char({ faults.name }),                        ...
             'FontSize', opt.FontSize, 'Color', opt.Color, ...
             'Interpreter', 'none');

   if nargout > 0, varargout{1} = {hf, ht}; end
end
