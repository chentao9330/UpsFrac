% test satnum



G = cartGrid([200, 1, 1]);
G = computeGeometry(G);

bc = addBC([], find(G.faces.centroids(:,1)==0), 'pressure', 1*barsa, 'sat', [1 0]);
bc = addBC(bc, find(G.faces.centroids(:,1)==G.cartDims(1)), 'pressure', 0, 'sat', [0 0]);

x = linspace(0, 1, 11);
y = linspace(1, 0, 11);

fluidSimple = initSimpleFluid('mu' , [   1,  10]*centi*poise     , ...
                              'rho', [1014, 859]*kilogram/meter^3, ...
                              'n'  , [   2,   2]);

kr = fluidSimple.relperm(x', initResSol(G, 0.0));

% linear relperm
grdecl_fluid.swof{1} = [x'  x' y' y'];      
grdecl_fluid.swof{2} = [x' kr(:,1) kr(:,2) y'];

satnum = ones(G.cells.num,1);
satnum(G.cells.num/2+1:end) = 2;


fluid  = initSatnumFluid(grdecl_fluid, 'mu' , [   1,  10]*centi*poise, ...
                              'rho', [1014, 859]*kilogram/meter^3, ...
                           'satnum', satnum, ...
                           'imbnum', satnum);
                        
%fluid = initSimpleFluid('mu', [1 1],  'n', [1 1]);
rSol = initResSol(G, 0, 0.0);

% add satnum fields to rSol
rSol.minSat = rSol.s;

rock.perm = ones(G.cells.num, 1);
rock.poro = 0.5*ones(G.cells.num, 1);
S = computeMimeticIP(G, rock);

%% Solve initial pressure
% Solve linear system construced from S and W to obtain solution for flow
% and pressure in the reservoir 
rSol = solveIncompFlow(rSol, G, S,  fluid, 'bc', bc);

%return
%% Main loop

T      = 1*milli*second(); %minute(); %day(); %year();
dT     = T/100;
pv     = poreVolume(G,rock);

sat_t = [];
kr_w = [];
kr_o = [];

% Prepare plotting of saturations
% clf
%    plotGrid(G,'FaceColor','none','EdgeAlpha',0.1);
%    
%    axis off, %view(30,50), 
%   % colormap(flipud(jet))
%    colorbar('horiz'); 
%    hs = []; ha=[];% zoom(2.5)
%    delete([hs, ha])
%    hs = plotCellData(G, rSol.s);
%    %[X,Y,Z]=meshgrid(1:46,1:112,3:4); plotCellData(G,rSol.s,cart2active(G,sub2ind(G.cartDims,X(:),Y(:),Z(:))))
%    %hs = plotCellData(G, rSol.s,find(rSol.s>0.01));
%    ha = annotation('textbox',[0.6 0.2 0.5 0.1], 'LineStyle','none', ...
%       'String', ['Water saturation at ',num2str(convertTo(0,year)),' years']);
%   axis tight off
  %drawnow
   
% Start the main loop
t  = 0;  plotNo = 1;
figure;

while t < T,
   rSol = explicitTransport(rSol, G, dT, rock, fluid, 'bc', bc, ...
                            'verbose', false, 'dt_factor', 0.1);

   % Check for inconsistent saturations
   assert(max(rSol.s(:,1)) < 1+eps && min(rSol.s(:,1)) > -eps);

   % Update solution of pressure equation.
   rSol = solveIncompFlow(rSol, G, S, fluid, 'bc', bc);

   % Increase time and continue if we do not want to plot saturations
   t = t + dT;
   if ( t < plotNo*dT  && t <T), continue, end

   %%
   % Plot saturation
%    figure(3)
%    delete([hs, ha])
%    hs = plotCellData(G, rSol.s,find(rSol.s>0.01));
%     % plotCellData(G, rSol.s);
%    ha = annotation('textbox',[0.6 0.2 0.5 0.1], 'LineStyle','none', ...
%      'String', ['Water saturation at ',num2str(convertTo(t,year)),' years']);
%    %view(30, 50+7*(plotNo-1)), 
%    drawnow
%    plotNo = plotNo+1;
%    
%    sat_t = [sat_t rSol.s];

   plot(rSol.s)
   ylim(gca, [ 0 1])
   ylabel('Saturation')
   xlabel('cell')
   title(['t = ', num2str(t), ' [seconds]']);
   drawnow;
end
