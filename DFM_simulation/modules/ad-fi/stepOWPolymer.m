function [state, meta] = stepOWPolymer(state0, state, meta, dt, G, W, system, fluid, varargin)
% Do a single step of a nonlinear solve for a Oil-Water-Polymer system
% This function should in general not be called directly and is as such not
% documented with regards to input/output: See solvefiADI for an
% explanation of how the ad-fi solvers are implemented.

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

% $Date$
% $Revision$

opt = struct('Verbose', mrstVerbose);
opt = merge_options(opt, varargin{:});
s = system.s;



if ~isfield(state, 'c')
    state.c = 0*state.pressure;
end

if ~isfield(state, 'cmax')
    state.cmax = state.c;
end

eqs = eqsfiOWPolymerExplicitWells(state0, state, dt, G, W, s, fluid);

if system.nonlinear.cpr && isempty(system.podbasis)
    [dx, gmresits, gmresflag] = cprGeneric(eqs, system,...
                                'ellipSolve', system.nonlinear.cprEllipticSolver,...
                                'cprType',    system.nonlinear.cprType,...
                                'relTol',     system.nonlinear.cprRelTol);
else
    dx = SolveEqsADI(eqs, system.podbasis);
end


[meta, residuals] = getResiduals(meta, eqs, system, 0);
[dx, meta] = stabilizeNewton(dx, meta, system);



searchfail = true;
if system.nonlinear.linesearch
    getEqs = @(state) eqsfiOWPolymerExplicitWells(state0, state, dt, G, W, s, fluid, 'resOnly', true);
    upState = @(dx) updateState(state, dx, fluid, W);
    [state, dx, searchfail] = linesearchADI(state, dx, system, getEqs, upState, false);
end

% Update reservoir conditions once a delta has been found.
if searchfail
    [dx, meta] = stabilizeNewton(dx, meta, system);
    % If the line search failed, uncritically accept the first step and
    % pray the other measures (relaxation / dampening) handle the error.
    [state, nInc] = updateState(state, dx, fluid, W);
end



[converged CNV MB] = getConvergence(state, eqs, fluid, system, dt);
meta.converged = converged;
meta.stopped = meta.iteration == system.nonlinear.maxIterations && ~converged;

if opt.Verbose
    residuals = cellfun(@(x) norm(x.val, 'inf'), eqs);
    eqnnames = {'Oil', 'Water', 'Polymer', 'qOs', 'qWs', 'pBHP'};
    printResidual(residuals, [], eqnnames, meta.iteration, CNV, MB);
end
end

%--------------------------------------------------------------------------

function [state, nInc] = updateState(state, dx, fluid, W)
maxSatStep = .2;



dp = dx{1};
ds = dx{2};
dc = dx{3};

nInc = max( norm(dp,'inf')/norm(state.pressure, 'inf'), ...
            norm(ds,'inf')/norm(state.s(:,1), 'inf') );
    
maxch = norm(ds, 'inf');
step = min(1, maxSatStep./maxch);


state.pressure = state.pressure + step*dp;
sw = state.s(:,1) + step*ds;

state.s = [sw, 1-sw];
state.c = state.c + dc;
state.c = min(state.c, fluid.cmax);
state.c = max(state.c, 0);


state.cmax = max(state.cmax, state.c);

dqWs  = dx{4};
dqOs  = dx{5};
% Intentionally skipping 6 - the sixth equation is a trivial addition which
% only serves a purpose for computations of the adjoint gradient. The
% result is already known in forward simulations.
dpBHP = dx{7};

for w = 1:numel(state.wellSol)
    state.wellSol(w).pressure = state.wellSol(w).pressure + dpBHP(w);
    state.wellSol(w).qWs      = state.wellSol(w).qWs + dqWs(w);
    state.wellSol(w).qOs      = state.wellSol(w).qOs + dqOs(w);
    if isfield(W, 'poly')
        poly = W(w).poly;
        if isempty(poly)
            state.wellSol(w).poly = 0;
        else
            state.wellSol(w).poly = W(w).poly;
        end
    end
end
end
