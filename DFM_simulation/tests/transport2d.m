% test satnum

clf
gravity off
verbose=mrstVerbose();
% Domain.nodes = [0,0;1,0;1,1;0,1];
% Domain.edges = [1,2;2,3;3,4;4,1];
%G = triangleGrid(Domain.nodes, Domain.edges, 'maxArea', 0.01);
G = cartGrid([50, 50, 1]);
G = computeGeometry(G);

rock.perm = repmat(100*milli*darcy, [G.cells.num, 1]);
rock.poro = repmat(0.5,             [G.cells.num, 1]);
pv        = poreVolume(G, rock);

src = addSource([], [1 G.cells.num], [1 -1]*sum(pv)/year(), ...
                'sat', [1 0; 0 0]);
bc  = pside([], G, 'LEFT', 10, 'sat', [1,0]); 

fluid = initSimpleFluid('mu' , [   1,  10]*centi*poise     , ...
                        'rho', [1014, 859]*kilogram/meter^3, ...
                        'n'  , [   4,   4]);
S = computeMimeticIP(G, rock);

xe = initResSol(G, 0, 0.0);
xe = solveIncompFlow(xe, G, S,  fluid, 'src',src, 'bc', bc);
xi = xe;

T  = 0.07*year();
dT = T/10;

% Start the main loop
t  = 0; 

figure(1)
subplot(1,3,1);
   hpe = plotCellData(G, convertTo(xe.pressure(1:G.cells.num), barsa));
   title('Pressure field [bar]');
   axis equal tight;

while t < T,
 
   t0 = tic;
   xe = explicitTransport(xe, G, dT, rock, fluid, 'src', src, 'bc', bc, ...
                         'verbose', verbose);
   toc(t0)
%%{
   t0 = tic;
   xi = implicitTransport(xi, G, dT, rock, fluid, 'src', src, 'bc', bc, ...
                          'verbose', verbose);
   toc(t0)
   %}
   t = t + dT;

   % Plot saturation
   figure(1)
   subplot(1,3,2)
      hse = plotCellData(G, xe.s(:,1));
      title('Explicit transport');
      axis equal tight

   subplot(1,3,3);
      hsi = plotCellData(G, xi.s(:,1));
      title('Implicit transport');
      axis equal tight;

   drawnow
end
