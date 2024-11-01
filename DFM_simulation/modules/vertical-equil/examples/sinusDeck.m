function deck = sinusDeck(cartdims, nsteps, dt, L, dy, H, theta, depth, phi, perm, rate, p_press)

% ;
% cartdims = opt.cartDims
nc=prod(cartdims);
    % define runspec
    deck.RUNSPEC.cartDims=cartdims;
    deck.RUNSPEC.DIMENS=cartdims;
    deck.RUNSPEC.OIL=1;
    %only for eclipse
    deck.RUNSPEC.WATER=1;
    %deck.RUNSPEC.FIELD=1;
    deck.RUNSPEC.METRIC=1; 
    deck.RUNSPEC.TABDIMS=[1     1    20    50    20    50     1    20    20     1    10     1    -1     0     1];
    deck.RUNSPEC.WELLDIMS=[5 10 2 1 5 10 5 4 3 0 1 1];
    deck.RUNSPEC.AQUDIMS=[0 0 0 0 10 10 0 0];
    deck.RUNSPEC.START=734139;
    % one sat num region
    deck.REGIONS.SATNUM=ones(nc,1);
    %define props    
    %deck.PROPS.SWOF{1}=[s,s.^alpha,(1-s).^alpha,s*0];
    % set up fluid but is not used in a vertical equilibrium calculation
    chop = @(x) min(max(0,x),1);
    s_wc = 0.0; %0.2;
    s_or = 0.0; %0.2;
    s = linspace(s_wc, 1 - s_or, 2)'; alpha = 2;
    s_star = (s - s_wc)/(1 - s_wc - s_or);
    swof = [s, chop(s_star.^alpha), chop((1-s_star).^alpha), s*0.0];
    %swof = [s, chop(s_star.^alpha), chop((1-s_star).^alpha), s*drho*norm(g)*0.0];
    %swof = [swof; [1.0 1.0 0.0 0.0]];
    deck.PROPS.SWOF{1} = swof;
    
    pres = convertTo(convertFrom(6000, psia), barsa);
    deck.SOLUTION.PRESSURE=ones(nc,1)*pres;
     
    %deck.PROPS.DENSITY=[900 1000 0.044000000000000];
    deck.PROPS.DENSITY = [600 1000 1];
    deck.PROPS.ROCK=[100 1.0e-6 NaN NaN NaN NaN];
    %deck.PROPS.ROCK=[4000 3.000000000000000e-06 NaN NaN NaN NaN];
    deck.PROPS.PVTW=[100 1.0 0.0 0.40 0];
    %deck.PROPS.PVDO{1} = [ [300, 800, 8000]'*barsa, [1.05, 1.02, 1.01]', [2.85, 2.99, 3]' ];
    deck.PROPS.PVDO{1}=[100 1 0.1;1000 1 0.1];
    % define summary
    deck.SUMMARY=[];
    % difine SC
    %%
    % define grid
    deck.SCHEDULE.step.val=ones(nsteps,1)*dt/day;
    %deck.SCHEDULE.step.control=ones(nsteps,1);
    deck.GRID=grdeclSloping(cartdims,[L dy H],'theta',theta,'amp',H/5,'lambda',L/4);
    deck.GRID.ZCORN=deck.GRID.ZCORN+depth;
    deck.GRID.ACTNUM=int32(ones(nc,1));
    
    deck.GRID.PORO=ones(nc,1)*phi;
    deck.GRID.PERMX=ones(nc,1)*perm;
    deck.GRID.PERMY=ones(nc,1)*perm;
    deck.GRID.PERMZ=ones(nc,1)*perm;

    
    %
   
    sr=0.3;sw=0.3;
    deck.SOLUTION.SWAT=ones(nc,1);
    deck.SOLUTION.SOIL=ones(nc,1).*0.0;
    % define needed quantitites for simulation
    mu=[deck.PROPS.PVDO{1}(1,3),deck.PROPS.PVTW(1,4)]*centi*poise();
    rho=deck.PROPS.DENSITY(1:2);
    deck.SCHEDULE.control.WELSPECS=...
    {...
    'I01'    'W'    [  ceil(cartdims(1)/2)]    [ ceil(cartdims(2)/2) ]    [1000]    'OIL'    [0]    'STD'    'SHUT'    'NO'    [0]    'SEG'    [0];...
    'P01'    'W'    [  cartdims(1)]    [cartdims(2)]                     [1004.03647]    'OIL'    [0]    'STD'    'SHUT'    'NO'    [0]    'SEG'    [0];...
    };%#ok
    radius=0.01;
    deck.SCHEDULE.control.COMPDAT=...
    {...
    'I01'     [  ceil(cartdims(1)/4)]    [ ceil(cartdims(2)/2) ]   [1]    [cartdims(3)]    'OPEN'    [0]    [0]    [radius]    [-1]    [0]    'Default'    'Z'    [-1];...
    'P01'    [  cartdims(1)]    [cartdims(2)]     [1]    [cartdims(3)]                                         'OPEN'    [0]    [0]    [radius]    [-1]    [0]    'Default'    'Z'    [-1];...
    };%#ok

% use scaled spe10 rates
deck.SCHEDULE.control.WCONINJE=...
    {...
    'I01'  'OIL'  'OPEN'  'RESV'  [rate]  [rate]  [500-89.4018]  [Inf]  [0]  [0]...
    };%#ok
deck.SCHEDULE.control.WCONPROD=...
    {...
    'P01'  'OPEN'  'BHP'  [Inf]  [Inf]  [Inf]  [Inf]  [Inf]  [p_press-89.4018]  [0]  [0]  [0];...
    };%#ok
deck.SCHEDULE.control=[deck.SCHEDULE.control;deck.SCHEDULE.control];
deck.SCHEDULE.control(2).WCONINJE{3}='SHUT';
% define grid
deck.SCHEDULE.step.val=ones(nsteps,1)*dt/day;
deck.SCHEDULE.step.val=[deck.SCHEDULE.step.val,deck.SCHEDULE.step.val*40];
