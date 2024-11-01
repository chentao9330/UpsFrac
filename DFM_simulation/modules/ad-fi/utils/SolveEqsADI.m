function dx = SolveEqsADI(eqs, phi)
useBasis = ~isempty(phi);

if useBasis
    % We have been provided a pod basis, use it to reduce equation sets
    % before solving.
    for i = 1:numel(eqs)
        eqs{i}.val = phi.basis{i}'*eqs{i}.val;
    end
end

numVars = cellfun(@numval, eqs)';
cumVars = cumsum(numVars); 
ii = [[1;cumVars(1:end-1)+1], cumVars];

eqs = cat(eqs{:});
tic
J = -eqs.jac{:};

if useBasis
    blkphi = blkdiag(phi.basis{:});
    J = blkphi'*J*blkphi;
end

tmp = J\eqs.val;

eqn = size(ii,1);
dx = cell(eqn,1);
for i = 1:eqn
    dx{i} = tmp(ii(i,1):ii(i,2));
end

if useBasis
    for i = 1:numel(dx)
        dx{i} = phi.basis{i}*dx{i};
    end
end
