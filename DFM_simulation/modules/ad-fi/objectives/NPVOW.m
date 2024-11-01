function obj = NPVOW(G, wellSols, schedule, varargin)
% Compute net present value of a schedule with well solutions

opt     = struct('OilPrice',             1.0 , ...
                 'WaterProductionCost',  0.1 , ...
                 'WaterInjectionCost',   0.1 , ...
                 'DiscountFactor',       0.0 , ...
                 'ComputePartials',      false, ...
                 'tStep' ,               []);
opt     = merge_options(opt, varargin{:});    
             
ro  = opt.OilPrice            / stb;
rw  = opt.WaterProductionCost / stb;
ri  = opt.WaterInjectionCost  / stb;
d   = opt.DiscountFactor;             


% pressure and saturaton vectors just used for place-holding
p  = zeros(G.cells.num, 1);
sW = zeros(G.cells.num, 1);

dts   = schedule.step.val;

if isempty(opt.tStep) %do all
    time = 0;
    numSteps = numel(dts);
else
    time = sum(dts(1:(opt.tStep-1)));
    numSteps = 1;
    dts = dts(opt.tStep);
end
    
obj = repmat({[]}, numSteps, 1);

for step = 1:numSteps
    sol = wellSols{step};
    nW  = numel(sol);
    pBHP = zeros(nW, 1); %place-holder
    qWs  = vertcat(sol.qWs);
    qOs  = vertcat(sol.qOs);
    
    if opt.ComputePartials
        [~, ~, qWs, qOs, ~] = initVariablesADI(p, sW, qWs, qOs, pBHP);
    end
    
    dt = dts(step);
    time = time + dt;
    
    injInx  = (vertcat(sol.sign) > 0);
    prodInx = ~injInx;
    obj{step} = ( dt*(1-d)^(time/year) )*...
                spones(ones(1, nW))*( (-ro*prodInx).*qOs ...
                             +(rw*prodInx - ri*injInx).*qWs );
end
