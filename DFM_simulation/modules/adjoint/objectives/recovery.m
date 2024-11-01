function [obj] = recovery(G, S, W, rock, fluid, simRes, schedule, controls, varargin)
% recovery -- objective function calculating water volume at last time step
%
% SYNOPSIS:
%   obj = (G, S, W, rock, fluid, simRes, schedule, controls, varargin)
%
% DESCRIPTION:
%   Computes value of objective function for given simulation, and partial
%   derivatives of variables if varargin > 6
% PARAMETERS:
%   simRes      -
%   
%
% RETURNS:
%   obj         - structure with fields
%        val    - value of objective function
%
%
%
%
% SEE ALSO:

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

numSteps = numel(simRes);
porv     = G.cells.volumes.*rock.poro;
obj.val  = sum( porv.*simRes(numSteps).resSol.s );

if nargin > 6
    partials = repmat( struct('v', [], 'p', [], 'pi', [], 'q_w', [], 's', [], 'u', []), [numSteps 1] );
    numCF    = size(G.cells.faces, 1);
    numC     = G.cells.num;
    numF     = G.faces.num;
    numU     = numel(controls.well);

    for k = 1 : numSteps
        ts = k;
        partials(ts).v    = zeros(1, numCF);
        partials(ts).p    = zeros(1, numC);
        partials(ts).pi   = zeros(1, numF);
        partials(ts).q_w  = zeros(1, length( vertcat(W.cells) ));
        if ts == numSteps
            partials(ts).s = porv';
        else
            partials(ts).s = zeros(1, numC);
        end
        partials(ts).u  = zeros(1, numU);
    end
    obj.partials = partials;
end
