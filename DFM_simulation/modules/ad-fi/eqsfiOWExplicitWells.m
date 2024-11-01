function [eqs, hst] = eqsfiOWExplicitWells(state0, state, dt, G, W, s, f, varargin)

opt = struct('Verbose', mrstVerbose, ...
             'reverseMode', false,...
             'scaling', [],...
             'resOnly', false,...
             'history', []);

opt = merge_options(opt, varargin{:});

if ~isempty(opt.scaling)
    scalFacs = opt.scaling;
else
    scalFacs.rate = 1; scalFacs.pressure = 1;
end

hst = opt.history;

% current variables: ------------------------------------------------------
p    = state.pressure;
sW   = state.s(:,1); 
pBHP = vertcat(state.wellSol.pressure);
qWs  = vertcat(state.wellSol.qWs);
qOs  = vertcat(state.wellSol.qOs);

% previous variables ------------------------------------------------------
p0  = state0.pressure;
sW0 = state0.s(:,1);
%--------------------------------------------------------------------------


%Initialization of independent variables ----------------------------------

zw = zeros(size(pBHP));
if opt.resOnly
    % ADI variables aren't needed since we are only computing the residual.
elseif ~opt.reverseMode
    [p, sW, qWs, qOs, pBHP]  = initVariablesADI(p, sW, qWs, qOs, pBHP);
else
    [p0, sW0, ~, ~, zw] = initVariablesADI(p0, sW0,...
                                            zeros(size(qWs)), ...
                                            zeros(size(qOs)), ...
                                            zeros(size(pBHP))...
                                            );
end

g  = norm(gravity);
[Tw, dzw, Rw, wc, perf2well, pInx, iInxW] = getWellStuff(W);

%--------------------
%check for p-dependent tran mult:
trMult = 1;
if isfield(f, 'tranMultR'), trMult = f.tranMultR(p); end

%check for p-dependent porv mult:
pvMult = 1; pvMult0 = 1;
if isfield(f, 'pvMultR') 
    pvMult =  f.pvMultR(p); 
    pvMult0 = f.pvMultR(p0);
end

%check for capillary pressure (p_cow)
pcOW = 0;
if isfield(f, 'pcOW') 
    pcOW  = f.pcOW(sW);
    pcOWw = pcOW(wc);
end

% -------------------------------------------------------------------------
% [krW, krO] = f.relPerm(sW);
krW = f.krW(sW);
krO = f.krOW(1-sW);

% water props (calculated at oil pressure OK?)
%bW     = f.bW(p);
bW     = f.bW(p-pcOW);
rhoW   = bW.*f.rhoWS;
% rhoW on face, avarge of neighboring cells (E100, not E300)
rhoWf  = s.faceAvg(rhoW);
%mobW   = trMult.*krW./f.muW(p);
mobW   = trMult.*krW./f.muW(p-pcOW);
dpW     = s.grad(p-pcOW) - g*(rhoWf.*s.grad(G.cells.centroids(:,3)));
% water upstream-index
upc = (double(dpW)>=0);
bWvW = s.faceUpstr(upc, bW.*mobW).*s.T.*dpW;


% oil props
bO     = f.bO(p);
rhoO   = bO.*f.rhoOS;
rhoOf  = s.faceAvg(rhoO);
mobO   = trMult.*krO./f.BOxmuO(p);
dpO    = s.grad(p) - g*(rhoOf.*s.grad(G.cells.centroids(:,3)));
% oil upstream-index
upc = (double(dpO)>=0);
bOvO   = s.faceUpstr(upc, bO.*mobO).*s.T.*dpO;


%WELLS ----------------------------------------------------------------
bWw     = bW(wc);
bOw     = bO(wc);
mobWw  = mobW(wc);
mobOw  = mobO(wc);

%producer mobility
bWmobWw  = bWw.*mobWw;
bOmobOw  = bOw.*mobOw;

%set water injector mobility: mobw = mobw+mobo+mobg, mobo = 0;
bWmobWw(iInxW) = bWw(iInxW).*(mobWw(iInxW) + mobOw(iInxW));
bOmobOw(iInxW) = 0;

pw  = p(wc);

bWqW  = -bWmobWw.*Tw.*(pBHP(perf2well) - pw + pcOWw + g*dzw.*rhoW(wc));
bOqO  = -bOmobOw.*Tw.*(pBHP(perf2well) - pw + g*dzw.*rhoO(wc));

% EQUATIONS ---------------------------------------------------------------


% oil:
eqs{1} = (s.pv/dt).*( pvMult.*bO.*(1-sW) - pvMult0.*f.bO(p0).*(1-sW0) ) + s.div(bOvO);
eqs{1}(wc) = eqs{1}(wc) + bOqO;

% water:
eqs{2} = (s.pv/dt).*( pvMult.*bW.*sW - pvMult0.*f.bW(p0).*sW0 ) + s.div(bWvW);
eqs{2}(wc) = eqs{2}(wc) + bWqW;

% well equations
zeroW = 0*zw;

eqs{3} = Rw'*bWqW + qWs + zeroW;
eqs{4} = Rw'*bOqO + qOs + zeroW;

% Last eq: boundary cond
eqs{5} = handleBC(W, pBHP, qWs, qOs, [], scalFacs) + zeroW;
end
%--------------------------------------------------------------------------










