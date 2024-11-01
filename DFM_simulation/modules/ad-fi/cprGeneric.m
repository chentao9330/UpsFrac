function [dx, its, fl] = cprGeneric(eqs, system, varargin)
% A generic CPR preconditioner for implicit systems.
% Currently assumes that
% - The first variable is the pressure
% - That any cellwise variables to be solved for are next.
% - That any variables which can be eliminated before solving comes after
% the cellwise variables.
% - Well closure equations are last.

opt = struct('cprType'    ,'colSum', ...
             'ellipSolve' ,@mldivide, ...
             'relTol'     , 2e-2,...
             'pscale'     , system.pscale);
         
opt = merge_options(opt, varargin{:});
ne = numel(eqs);

copt = system.cpr;
active = max(copt.active);
g = copt.gas;
if ~isempty(g)
    % for unsaturated cells, switch respective columns of dx/ds and dx/drs
    unSat = logical(full(diag(eqs{g(2)}.jac{g(1)})));
    eqs   = switchCols(eqs, g, unSat);
end

n_elim = ne-active;
eliminated = cell(n_elim,1);

for i = 1:n_elim
    [eqs, eq] = elimVars(eqs, active + 1);
    eliminated{i} = eq;
end

ii = getEquationInxs(eqs);

if system.nonlinear.cprBlockInvert
    [A, b, pInx] = getCPRSystemBlockInvert(eqs, ii, opt, active);
else
    [A, b, pInx] = getCPRSystemDiagonal(eqs, ii, opt);
end

% Scale pressure variables
pscale = opt.pscale;
if pscale ~= 1;
    pind = ii(1,1):ii(1,2);
    A(:,pind) = A(:,pind) ./ pscale;
end

prec           = getTwoStagePrec(A, pInx, opt);


[dX, fl, relres, its] = gmres(A, b, [], opt.relTol, 20, prec);

if pscale ~= 1;
    dX(pind) = dX(pind) ./ pscale;
end

dx = cell(ne, 1);
for i = 1:active
    dx{i} = dX(ii(i,1):ii(i,2));
end

% Recover eliminated variables. This is done in the reverse order of the
% elimination.

for i = n_elim:-1:1
    dVal = recoverVars(eliminated{i}, active + 1, {dx{[1:active (active + i + 1 : ne)]}}); %#ok<CCAT1>
    dx{active + i} = dVal;
end


%Assign dsG and dRS
if ~isempty(g)
    dM1 = dx{g(1)};
    dM2 = dx{g(2)};
    dx{g(1)} = ~unSat.*dM1  +  unSat.*dM2;
    dx{g(2)} =  unSat.*dM1  + ~unSat.*dM2;
end

if fl ~= 0
    warning('GMRES did not converge, Relative residual: %9.2e, error code: %2d', relres, fl);
end
end

function eInx = getEquationInxs(eqs)
numVars = cellfun(@numval, eqs)';
cumVars = cumsum(numVars); 
eInx = [[1;cumVars(1:end-1)+1], cumVars];
end

%--------------------------------------------------------------------------
function eqs   = switchCols(eqs, n, inx)
for k = 1:numel(eqs)
    tmp = eqs{k}.jac{n(1)}(:, inx);
    eqs{k}.jac{n(1)}(:, inx) = eqs{k}.jac{n(2)}(:, inx);
    eqs{k}.jac{n(2)}(:, inx) = tmp;
end 
end

%--------------------------------------------------------------------------
function [eqs, eqn] = elimVars(eqs, n)
% eliminate set of unknowns nr n using equation n () 
solveInx = setdiff(1:numel(eqs), n);
eqn      = eqs{n};

for eqNum = solveInx
    for jacNum = solveInx
        eqs{eqNum}.jac{jacNum} = eqs{eqNum}.jac{jacNum} - eqs{eqNum}.jac{n}*(eqn.jac{n}\eqn.jac{jacNum});
    end
    eqs{eqNum}.val = eqs{eqNum}.val - eqs{eqNum}.jac{n}*(eqn.jac{n}\eqn.val);
end

eqs  = eqs(solveInx);
for eqNum = 1:numel(eqs)
    eqs{eqNum}.jac = eqs{eqNum}.jac(solveInx);
end

end
%--------------------------------------------------------------------------
function x = recoverVars(eq, n, sol)


% recover variables x at position n using solutions sol
solInx = [1:(n-1), (n+1):(numel(sol)+1)];
x = - eq.jac{n}\(eq.val);
for k  = 1:numel(solInx)
    x = x - eq.jac{n}\(eq.jac{solInx(k)}*sol{k});
end
end

function [A, b, pInx] = getCPRSystemDiagonal(eqs, ii, opt)
pInx = false(ii(end,end), 1);
pInx(ii(1,1):ii(1,2)) = true;


if strcmpi(opt.cprType, 'diag')
    cprFunc = @diag;
elseif strcmpi(opt.cprType, 'colsum')
    cprFunc = @sum;
end


n = numel(eqs{1}.val);

deqs = eqs;
for k = 1:numel(eqs)
    for l = 1:numel(eqs{1}.jac)
        deqs{k}.jac{l} = spdiags(cprFunc(eqs{k}.jac{l})', 0, n,n);
    end
end
deqs = cat(deqs{:});
D = deqs.jac{:};

eqs = cat(eqs{:});
A   = eqs.jac{:};
b   = -eqs.val;


A = D\A;
b = D\b;
end

function [A, b, pInx] = getCPRSystemBlockInvert(eqs, ii, opt, active)
assert(active == 2 || active == 3);

pInx = false(ii(end,end), 1);
pInx(ii(1,1):ii(1,2)) = true;



if strcmp(opt.cprType, 'diag')
    cprFunc = @diag;
elseif strcmp(opt.cprType, 'colSum')
    cprFunc = @sum;
end

if active == 3
    % Solving for pressure and two other cell wise variables
    dss = cell(2,2);
    dps = cell(1,2);
    pI = 1;
    sI = [2,3];
    for ii = 1:2
        for jj = 1:2
            dss{ii,jj} = reshape(cprFunc(eqs{sI(ii)}.jac{sI(jj)}), [], 1); 
        end
        dps{ii} = reshape(cprFunc(eqs{pI}.jac{sI(ii)}), [], 1); 
    end

    n = numel(dss{1,1});
    inx = (1:n)';
    dtrmInv = 1./(dss{1,1}.*dss{2,2} - dss{2,1}.*dss{1,2});

    DssInv  = sparse([inx, inx  , inx+n, inx+n], ...
                     [inx, inx+n, inx  , inx+n], ...
                     [dss{2,2}, -dss{1,2}, -dss{2,1}, dss{1,1}].*(dtrmInv*[1 1 1 1]), ...
                     2*n, 2*n);


    Dps     = sparse([inx, inx], [inx, inx+n], [dps{1}, dps{2}], n, 2*n);
else
    % A simpler two variable system
    pI = 1;
    sI = 2;
    dss = reshape(cprFunc(eqs{sI}.jac{sI}), [], 1); 
    dps = reshape(cprFunc(eqs{pI}.jac{sI}), [], 1); 

    n = numel(dss);
    inx = (1:n)';

    DssInv  = sparse(inx, ...
                     inx, ...
                     1./dss, ...
                     n,n);

    Dps     = sparse(inx, inx, dps, n, n);
end


M = Dps*DssInv;

eqs = cat(eqs{:});
A   = eqs.jac{:};
b   = -eqs.val;

A(pInx, :) = A(pInx, :) - M*A(~pInx,:);
b(pInx)    = b(pInx)    - M*b(~pInx);
end

%--------------------------------------------------------------------------

function prec = getTwoStagePrec(A, pInx, opt)
Ap     = A(pInx, pInx);

setup.type = 'nofill';
[L, U] = ilu(A, setup);

prec   = @(r)applyTwoStagePrec(r, A, L, U, Ap, pInx, opt);
end

%--------------------------------------------------------------------------

function x = applyTwoStagePrec(r, A, L, U, Ap, pInx, opt)
x = zeros(size(r));
x(pInx) = opt.ellipSolve(Ap, r(pInx));

r = r - A*x;
x = x + U\(L\r);
end



