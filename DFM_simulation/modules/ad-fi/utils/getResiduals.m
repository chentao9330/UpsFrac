function [meta, residuals] = getResiduals(meta, eqs, system, gmresflag)
    % Store the residuals for debugging and convergence testing.
    residuals = cellfun(@(x) norm(x.val, 'inf'), eqs);

    if ~isfield(meta, 'res_history')
        meta.res_history = zeros(system.nonlinear.maxIterations, numel(residuals));
    end

    meta.res_history(meta.iteration, :) = residuals;

    % Try a simple detection of oscillations, and relax the next iteration if
    % oscillations were detected.
    if detectNewtonOscillations(meta.res_history, system.cellwise, meta.iteration, system.nonlinear.relaxRelTol) && ~gmresflag
        meta.relax = max(meta.relax - system.nonlinear.relaxInc, system.nonlinear.relaxMax);
        fprintf(['Oscillating behavior detected: Relaxation set to ' num2str(meta.relax) '\n'])
    end
end
