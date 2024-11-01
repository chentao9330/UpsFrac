% Files
%  SolveEqsADI              - We have been provided a pod basis, use it to reduce equation sets
%  computeNumGrad           - Compute numerical gradient w.r.t wells
%  dinterpReg               - if reginx{k} == ':' for improved eff seperate this case
%  getConvergence           - Compute convergence based on total mass balance (tol_mb) and maximum
%  getRegMap                - if nargin < 4 whole domain
%  getResiduals             - Store the residuals for debugging and convergence testing.
%  getSimMatrices           - half-trans -> trans and reduce to interior
%  getWellStuff             - ------------------------------------------
%  handleRegions            - also need actnum/G.cells.indexmap
%  initWellSolLocal         - Initialize well solution data structure.
%  interpReg                - if reginx{k} == ':' for improved eff seperate this case
%  printResidual            - fprintf('-9s', eqnnames{:})
%  scheduleFromSummary      - There are n-1 control steps for n summary steps
%  setupSimComp             - Set up helper structure for solvers based on automatic differentiation.
%  solveAdjointEqsADI       - If adjVec is not empty, we are not at the last timestep (first in the
%  updateSchedule           - --------------------------------------------------------------------------
%  wellSolToVector          - Helper function which makes cell arrays of well solutions easier to plot

%{
#COPYRIGHT#
%}
