close all
%clear all

%% Multiscale Pressure Solver with capillary pressure: 
%  Flow Driven by two Vertical Wells
%  
%
% 2D case with homogeneous permeability and porosity but 2 satnum
% regions where the capillary pressure curve is different in the two
% regions.

simpleGrid = false;
verbose    = true;

pc_form = 'nonwetting';

%% construct simple cartesian testcase
nx = 20; ny = 20; nz = 1;
G         = cartGrid([nx ny nz]);
G         = computeGeometry(G);
rock.perm = repmat(100*milli*darcy, [G.cells.num, 1]);
rock.poro  = repmat(0.3, [G.cells.num, 1]);

% define two satnum regions
satnum = ones(G.cells.num,1);
satnum((end/2)+1:end) = 2;

x = linspace(0, 1, 11)';
y = linspace(1, 0, 11)';

%% Make fluid tables (normally result of readpvt)
% define different capillary pressure in the different regions

props.swof{1} = [x x y y*0.2*barsa];  % s kr_w(s) kr_o(s) pc(s) 
props.swof{2} = [x x y  y*0.1*barsa]; 

% initialize fluid
fluid = initSWOFFluid('mu' , [   1,  10] .* centi*poise     , ...
                      'rho', [1000, 700] .* kilogram/meter^3, ...
                      'table', props.swof, ...
                      'satnum', satnum);
           

%% Set wells
rate = 0.5*meter^3/day;
bhp  = 1*barsa;

W = verticalWell([], G, rock, 1, 1, 1:nz,          ...
                  'Type', 'rate', 'Val', rate, ...
                  'Radius', .1, 'Name', 'I', 'Comp_i', [1 0]);
W = verticalWell(W, G, rock, nx, ny, 1:nz,     ...
                  'Type','bhp', 'Val', bhp, ...
                  'Radius', .1, 'Dir', 'x', 'Name', 'P', 'Comp_i', [0 1]);


%% Set up solution structures
% Here we need four solution structures, two for each simulator to hold the
% solutions on the grid and in the wells, respectively.

xRef = initState(G, W, 0, [0.2, 0.8]);

%% Assemble linear systems
gravity off
S  = computeMimeticIP(G, rock, 'Verbose', verbose);


%% Do some plotting of pc-curves
xDummy  = initResSol (G, 0.0, 0);
half = G.cells.num/2;
% to plot for different rock types, numel(s) must be g.cells.num
xDummy.s = [linspace(0, 1, numel(xDummy.s)/2)'; ...
   linspace(0, 1, numel(xDummy.s)/2)'];
pc = convertTo(fluid.pc(xDummy), barsa);

figure
subplot(1, 2, 1)
plot(xDummy.s(1:half), pc(1:half)); hold on
plot(xDummy.s(half+1:end), pc(half+1:end), 'r');
xlabel('s_w'); ylabel('pc [bar]'); legend('region 1', 'region 2');
subplot(1, 2, 2)
title('Regions'); plotCellData(G, satnum);
plotGrid(G, 'faceColor', 'none');


%% Set up pressure and transport solvers
psolve  = @(state) solveIncompFlow  (state, G, S, fluid, 'wells', W);
tsolve = @(state, dT) explicitTransport(state, G, dT, rock, fluid, 'wells', W, 'verbose', verbose); 
%tsolve = @(state, dT) implicitTransport(state, G, dT, rock, fluid, 'wells', W, 'verbose', verbose);


%% Solve initial pressure
xRef = psolve(xRef);

%% Plot initial pressure
figure
   plotCellData(G, convertTo(xRef.pressure(1:G.cells.num), barsa));
   title('Initial pressure fine'),
   cx = caxis;
   axis off; view(2)   

% Report pressure in wells.
dp = @(x) num2str(convertTo(x.wellSol(1).pressure-x.wellSol(2).pressure, barsa));
disp(['DeltaP, Fine: ', dp(xRef)])


%% Transport loop
   T      = 1*year;
   dT     = year/10;
   dTplot = dT*2;

N      = fix(T/dTplot);
pv     = poreVolume(G,rock);
%%
verbose = true;
tic

%% Start the main loop
t  = 0; plotNo = 1; hi = ''; he = '';
e = []; pi = []; pe = [];
figure;

while t < T,   
   % TRANSPORT SOLVE
   xRef = tsolve(xRef, dT);  
                             
   % Check for inconsistent saturations
   s = xRef.s(:,1);
   assert(max(s) < 1+eps && min(s) > -eps); 
   assert(~any(isnan(s)) && ~any(isinf(s))); 
      
   % Update solution of pressure equation.
   xRef =  psolve(xRef) ; %   
   
   % Measure water saturation in production cells in saturation  
   pe = [pe; xRef.s(W(2).cells,1)' ];                 %#ok
  
   % Increase time and continue if we do not want to plot saturations
   t = t + dT;   
   
   if ( t < plotNo*dTplot && t <T), continue, end
   % Plot saturation
   %%{
   r = 0.01;
   clf
   title(['Saturation at ', num2str(t/day), ' days']); 
   plotGrid(G, 'faceColor', 'none', 'edgeAlpha', 0.1);
   plotCellData(G, xRef.s(:,1), find(xRef.s>0.01)); caxis([0 1])
   view(2), 
   axis off

   drawnow
   plotNo = plotNo+1;
   %}
end


n = size(pe,1);
figure
   plot(1:n,pe(:,1),'-o')
   title('Water saturation in production well');


