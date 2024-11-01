% horizonalWellOverlap;
% Long horizontal well using block- and well overlap

% Clear
cellDims = [60, 60, 5];
resDims  = [500 500 20];
verbose  = true;
sizLayer = prod(cellDims(1:2));

blockOverlap = 6;
wellOverlap  = 6;

G = computeGeometry(cartGrid(cellDims, resDims));

% layered perm-field
lo = ones([sizLayer, 1]);
lp = exp(5*( rand(cellDims(3),1) - .5) );
rock.perm = convertFrom(kron(lp, lo), milli*darcy);

% 3 vert inj, 1 hor prod
W = addWell([], G, rock, 49*60+55 : sizLayer : G.cells.num,      ...
            'Type', 'rate', 'Val', 1/3*meter^3/day, 'Radius', 1, ...
            'InnerProduct', 'ip_tpf', 'Comp_i', [1, 0]);
W = addWell(W, G, rock, 1         : sizLayer : G.cells.num,      ...
            'Type', 'rate', 'Val', 1/3*meter^3/day, 'Radius', 1, ...
            'InnerProduct', 'ip_tpf', 'Comp_i', [1, 0]);
W = addWell(W, G, rock, 30        : sizLayer : G.cells.num,      ...
            'Type', 'rate', 'Val', 1/3*meter^3/day, 'Radius', 1, ...
            'InnerProduct', 'ip_tpf', 'Comp_i', [1, 0]);

W = addWell(W, G, rock, 2*sizLayer + 14*60 + (10:55),         ...
            'Type', 'bhp', 'Val', 0, 'Radius', 1, 'Dir', 'x', ...
            'InnerProduct', 'ip_tpf', 'Comp_i', [0, 1]);

part  = partitionUI(G, [4, 4, 1]);
part2 = processPartition  (G, part,  'Verbose', verbose);
CG    = generateCoarseGrid(G, part2, 'Verbose', verbose);

S = computeMimeticIP(G, rock,            ...
                     'Verbose', verbose, ...
                     'Type', 'comp_hybrid',   ...
                     'InnerProduct', 'ip_tpf');

CS = generateCoarseSystem (G, rock, S, CG, ones([G.cells.num, 1]), ...
                           'Verbose', verbose, 'Overlap', blockOverlap);

W = generateCoarseWellSystem(G, S, CG, CS, ones([G.cells.num, 1]), ...
                             rock, W, 'OverlapWell', wellOverlap,  ...
                             'OverlapBlock', blockOverlap);

stateRef = initState(G, W, 0, [0, 1]);
stateMs  = stateRef;

fluid    = initSimpleFluid('mu' , [   1,  10]*centi*poise     , ...
                           'rho', [1000, 700]*kilogram/meter^3, ...
                           'n'  , [   2,   2]);

stateRef = solveIncompFlow  (stateRef, G, S, fluid, ...
                             'wells', W, 'Solver', 'mixed');
stateMs  = solveIncompFlowMS(stateMs, G, CG, part2, S, CS, ...
                             fluid, 'wells', W, 'Solver', 'mixed');

disp(['DeltaP - Fine: ', num2str(convertTo(stateRef.wellSol(1).pressure, barsa))])
disp(['DeltaP - Ms:   ', num2str(convertTo(stateMs .wellSol(1).pressure, barsa))])

%% plot output
f = figure;
subplot(2,2,1)
   plotGrid(G, 'FaceColor', 'none', 'EdgeColor', [.8, .8, .8])
   plotGrid(G, reshape([W(1:3).cells], [], 1), 'FaceColor', 'b', ...
            'EdgeColor', 'black');
   plotGrid(G, W(4).cells, 'FaceColor', 'r', 'EdgeColor', 'black');
   axis tight
   title('Well-cells')

subplot(2,2,2)
   plot(-convertTo(stateRef.wellSol(4).flux, meter^3/day), 'b'); hold on
   plot(-convertTo(stateMs .wellSol(4).flux, meter^3/day), 'r');
   axis([0 46 0 .1]); 
   legend('Fine', 'Multiscale')
   title('Producer inflow profile')

subplot(2,2,3)
   plotCellData(G, convertTo(stateRef.pressure(1:G.cells.num), barsa));
   title('Pressure Fine [bar]')
   axis tight , cax = caxis; colorbar

subplot(2,2,4)
   plotCellData(G, convertTo(stateMs .pressure(1:G.cells.num), barsa));
   title('Pressure Coarse [bar]')
   axis tight, caxis(cax); colorbar
