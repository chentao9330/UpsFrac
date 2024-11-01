function D = computeTOFandTracer(state, G, rock,  varargin)
%Compute time-of-flight and tracer distribution using finite-volume scheme.
%
% SYNOPSIS:
%    D  = computeTOFandTracer(state, G, rock)
%    D  = computeTOFandTracer(state, G, rock, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   Construct the basis for flow diagnostic by computing
%     1) time-of-flight        :   \nabla·(v T) = \phi,
%     2) reverse time-of-flight:  -\nabla·(v T) = \phi,
%     3) stationary tracer     :  -\nabla·(v C) = 0
%   using a first-order finite-volume method with upwind flux.
%
% REQUIRED PARAMETERS:
%   G     - Grid structure.
%
%   rock  - Rock data structure.
%           Must contain a valid porosity field, 'rock.poro'.
%
%   state - Reservoir and well solution structure either properly
%           initialized from functions 'initResSol' and 'initWellSol'
%           respectively, or the results from a call to function
%           'solveIncompFlow'.  Must contain valid cell interface fluxes,
%           'state.flux'.
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%   wells - Well structure as defined by function 'addWell'.  May be empty
%           (i.e., wells = []) which is interpreted as a model without any
%           wells.
%
%   src   - Explicit source contributions as defined by function
%           'addSource'.  May be empty (i.e., src = []) which is
%           interpreted as a reservoir model without explicit sources.
%
%   bc    - Boundary condition structure as defined by function 'addBC'.
%           This structure accounts for all external boundary conditions
%           to the reservoir flow.  May be empty (i.e., bc = []) which is
%           interpreted as all external no-flow (homogeneous Neumann)
%           conditions.
% RETURNS:
%   D - struct that contains the basis for computing flow diagnostics:
%       'inj'     - list of injection wells
%       'prod'    - list of production wells
%       'tof'     - time-of-flight and reverse time-of-flight returned
%                   as an (G.cells.num x 2) vector
%       'itracer' - steady-state tracer distribution for injectors
%       'ipart'   - tracer partition for injectors
%       'ptracer' - steady-state tracer distribution for producers
%       'ppart'   - tracer partition for producers

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

% Process optional parameters
opt = struct('bc', [], 'src', [], 'wells', []);
opt = merge_options(opt, varargin{:});
assert(~isempty(opt.wells));
assert(isempty(opt.src), 'Source terms not supported yet');
assert(isempty(opt.bc),  'Boundary conditions not supported yet');
%assert (isempty(opt.src) || ~isempty(opt.src.sat), ...
%   'Source terms must have a valid ''sat'' field.');
%assert (isempty(opt.bc) || ~isempty(opt.bc.sat), ...
%   'Boundary conditions must have a valid ''sat'' field.');

assert (isfield(rock, 'poro')         && ...
        numel(rock.poro)==G.cells.num,   ...
        ['The rock input must have a field poro with porosity ',...
         'for each cell in the grid.']);
assert(min(rock.poro) > 0, 'Rock porosities must be positive numbers.');


% Find injectors and producers
iwells = cellfun(@sum, {state.wellSol.flux})>0;
D.inj  = find( iwells);
D.prod = find(~iwells);

% Compute time-of-flight and tracer partition from injectors
t = computeTimeOfFlight(state, G, rock, 'wells', opt.wells, ...
   'tracer', {opt.wells(D.inj).cells});
D.tof     = t(:,1);
D.itracer = t(:,2:end);
[~,I]     = sort(D.itracer,2,'descend');
D.ipart   = I(:,1);

% Compute time-of-flight and tracer partition from producers
t = computeTimeOfFlight(state, G, rock, 'wells', opt.wells, ...
   'tracer', {opt.wells(D.prod).cells}, 'reverse', true);
D.tof(:,2) = t(:,1);
D.ptracer  = t(:,2:end);
[~,I]      = sort(D.ptracer,2,'descend');
D.ppart    = I(:,1);
