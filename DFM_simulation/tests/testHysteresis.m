% Test hysteresis by first injecting gas into a water filled reservoir (1D)
% and then injecting water at t = T/2
% 
close all
n = 300;
verbose=mrstVerbose();
% Make grid, initialize solution structure and set rock data
G = cartGrid([n, 1, 1]); %, [2, 1, 1]);
G = computeGeometry(G);
rSol = initResSol(G, 0, 1);
rock.perm = ones(G.cells.num, 1)*milli*darcy;
rock.poro = 0.5*ones(G.cells.num, 1);

%% Define 'fake' grdecl
% i.e instead of reading a grid file, set variables explicitly:
grdecl.ACTNUM = ones(n,1);
grdecl.SATNUM = ones(n,1);
grdecl.IMBNUM = 2*ones(n,1);

% add satnumfields to rSol struct
rSol.extSat = [rSol.s(:,1),rSol.s(:,1)];           


%% Read relperm curves for drainage and imbibition and construct fluid
SWFN{1} = [...
    0.2000         0    0.000
    0.3000    0.0278    0.000
    0.4000    0.1111    0.000
    0.5000    0.2500    0.000
    0.6000    0.4444    0.000
    0.7000    0.6944    0.000
    0.8000    1.0000    0.000
    0.9000    1.0       0.000
    1.0000    1.0       0.000 ];
SWFN{2} = SWFN{1};
 
 SGFN{1} = [...
      0         0         0
    0.1000    0.0156    0.000
    0.2000    0.0625    0.000
    0.3000    0.1406    0.000
    0.4000    0.2500    0.000
    0.5000    0.3906    0.000
    0.6000    0.5625    0.000
    0.7000    0.7656    0.000
    0.8000    1.0000    0.000
    0.9000    1.0       0.000
    1.0000    1.0       0.000 ];
 SGFN{2} = [...
         0         0         0
    0.1000         0    0.000
    0.2000         0    0.000
    0.3000         0    0.000
    0.4000         0    0.000
    0.5000    0.0625    0.000
    0.6000    0.2500    0.000
    0.7000    0.5625    0.000
    0.8000    1.0000    0.000
    0.9000    1.0       0.000
    1.0000    1.0       0.000 ];
 

grdecl_fluid.swfn = SWFN; 
grdecl_fluid.sgfn = SGFN; 

fluid = initSatnumFluid(grdecl_fluid, ...
                        'mu', [ 1.0, 0.318].*centi*poise, ...
                        'rho', [1033, 859].*kilogram/meter^3, ...
                        'satnum', grdecl.SATNUM, 'imbnum', grdecl.IMBNUM );

%% Plot relperm curves and different scanning curves:
r = rSol; r.s =linspace(0, 1, n)';
figure;
r.extSat = [r.s,r.s]; 
kr =  fluid.relperm(r.s, r);
plot(r.s, kr);
hold on;

% imbibition curve:
r.extSat(:,1) = 0.2; kr = fluid.relperm(r.s, r); plot(r.s, kr(:,2), 'r'); 
r.extSat(:,1) = 0.3; kr = fluid.relperm(r.s, r); plot(r.s, kr(:,2), 'c'); 
% scanning curves:
r.extSat(:,1) = 0.5; kr = fluid.relperm(r.s, r); plot(r.s, kr(:,2), 'k'); 
legend('drainage water', 'drainage gas', 'imbibtion gas', ...
       'scan_{gas}: 0.3', 'scan_{gas}: 0.5');


%% Set bc and compute mimetic ip
bc = addBC([], find(G.faces.centroids(:,1)==0), 'pressure', 300*barsa, ...
          'sat', [0 1]);
bc = addBC(bc, find(G.faces.centroids(:,1)==G.cartDims(1)), 'pressure', ...
           0, 'sat', [1 0]);
S = computeMimeticIP(G, rock);

%% Solve initial pressure
% Solve linear system construced from S and W to obtain solution for flow
% and pressure in the reservoir 
rSol = solveIncompFlow(rSol, G, S, fluid, 'bc', bc);

%return
%% Main loop
T      = 1*year();
dT     = 50*day;
pv     = poreVolume(G,rock);

figure;
hold off;
% Start the main loop
t  = 0;  plotNo = 1;
%figure;
while t < T,
   rSol = implicitTransport(rSol, G, dT, rock, fluid, 'bc', bc, ...
                            'verbose', verbose); %, 'dt', 1*day, 'ComputeDt', false);
   %rSol = explicitTransport(rSol, G, dT, rock, fluid, 'bc', bc, ...
   %                         'verbose', verbose);
   % Check for inconsistent saturations
   assert(max(rSol.s(:,1)) < 1+eps && min(rSol.s(:,1)) > -eps);

   % Update solution of pressure equation.
   rSol = solveIncompFlow(rSol, G, S, fluid, 'bc', bc);

   % Change to water injection at t = T/2
   if t > 0.5*T,  bc(1).sat = [1 0];end

   plot(rSol.s)
   ylim(gca, [ 0 1])
   ylabel('Saturation')
   xlabel('cell')
   title(['t = ', num2str(convertTo(t, day)), ' [days]']);
   drawnow

   t = t + dT;
end
