% test satnum


close all
G = cartGrid([40, 40, 1]);
G = computeGeometry(G);

src = addSource([], [1 G.cells.num], [1 -1], 'sat', [1 0; 0 0]);

x = linspace(0, 1, 11);
y = linspace(1, 0, 11);

fluidSimple = initSimpleFluid('mu' , [   1,  10]*centi*poise     , ...
                              'rho', [1014, 859]*kilogram/meter^3, ...
                              'n'  , [   4,   4]);

kr = fluidSimple.relperm(x', initResSol(G, 0.0));

%% Make fluid tables (normally result of readpvt)
% linear relperm
grdecl_fluid.swof{1} = [x'  x' y' zeros(numel(y),1)];  
% quadratic relperm
grdecl_fluid.swof{2} = [x' kr(:,1) kr(:,2) zeros(numel(y),1)];

sat1 = (G.cells.centroids(:,2)<=G.cartDims(2)-G.cells.centroids(:,1));
satnum = ones(G.cells.num,1);
satnum(~sat1) = 2;


fluid  = initSatnumFluid(grdecl_fluid, 'mu' , [   1,  10]*centi*poise, ...
                              'rho', [1014, 859]*kilogram/meter^3, ...
                              'satnum', satnum, ...
                              'imbnum', satnum);

rSol = initResSol(G, 0, 0.0);
rSol.minSat = rSol.s;

% plot saturation at x = y
plotCells = G.cells.centroids(:,2)==G.cells.centroids(:,1);

rock.perm = ones(G.cells.num, 1);
rock.poro = 0.5*ones(G.cells.num, 1);
S = computeMimeticIP(G, rock);

%% Solve initial pressure
% Solve linear system construced from S and W to obtain solution for flow
% and pressure in the reservoir
rSol = solveIncompFlow(rSol, G, S,  fluid, 'src', src);

%return
%% Main loop
T      = 10*minute(); %day(); %year();
dT     = T/100;
pv     = poreVolume(G,rock);

% Prepare plotting of saturations
figure(1)
clf
   plotGrid(G,'FaceColor','none','EdgeAlpha',0.1);
   axis off, 
   % colormap(flipud(jet))
   colorbar('horiz'); 
   hs = []; ha=[];% zoom(2.5)
   delete([hs, ha])
   hs = plotCellData(G, rSol.s);
   ha = annotation('textbox',[0.6 0.2 0.5 0.1], 'LineStyle','none', ...
      'String', ['Water saturation at ',num2str(convertTo(0,year)),' years']);
  axis tight off
  drawnow
   
% Start the main loop
t  = 0;  plotNo = 1;

while t < T,
   rSol = explicitTransport(rSol, G, dT, rock, fluid, 'src', src, ...
                            'verbose', false, 'dt_factor', 0.1);
  
   % Check for inconsistent saturations
   assert(max(rSol.s(:,1)) < 1+eps && min(rSol.s(:,1)) > -eps);

   % Update solution of pressure equation.
   rSol = solveIncompFlow(rSol, G, S, fluid, 'src',src);

   % Increase time and continue if we do not want to plot saturations
   t = t + dT;
   if ( t < plotNo*dT  && t <T), continue, end

   %%
   % Plot saturation
   figure(1)
   delete([hs, ha])
   hs = plotCellData(G, rSol.s(:,1));
   ha = annotation('textbox',[0.6 0.2 0.5 0.1], 'LineStyle','none', ...
     'String', ['Water saturation at ',num2str(convertTo(t,year)),' years']);
   drawnow
   plotNo = plotNo+1;
       
   figure(2)
   plot(rSol.s(plotCells,1))
   ylim(gca, [ 0 1])
   title('Saturation at cells x=y')
end
