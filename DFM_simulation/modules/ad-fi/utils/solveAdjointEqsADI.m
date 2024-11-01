function [x, ii] = solveAdjointEqsADI(eqs, eqs_p, adjVec, objk)
numVars = cellfun(@numval, eqs)';
cumVars = cumsum(numVars); 
ii = [[1;cumVars(1:end-1)+1], cumVars];

if iscell(objk), objk = objk{:};end
objk = cat(objk);
rhs  = -objk.jac{:}';
if ~isempty(adjVec)
    % If adjVec is not empty, we are not at the last timestep (first in the
    % adjoint recurrence formulation). This means that we subtract
    % the previous jacobian times the previous lagrange multiplier.
    % Previous here means at time t + 1.
    eqs_p = cat(eqs_p{:});
    rhs = rhs - eqs_p.jac{:}'*adjVec;
end
eqs = cat(eqs{:});
tic
x = eqs.jac{:}'\rhs;
tt = toc;
dispif(false, 'Lin. eq: %6.5f seconds, ', tt);
