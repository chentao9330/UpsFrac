G = processGRDECL(simpleGrdecl([20,20,15], 0.1)); 
%disp('splitting faces');tic, G = splitFaces(G, 1:G.faces.num);toc

gravity off
G = computeGeometry(G);
rock.perm = ones(G.cells.num, 1)*100*milli*darcy();
rock.poro = 0.5*ones(G.cells.num, 1);
pv     = poreVolume(G,rock);

src   = [];
%src   = addSource(src, [1 G.cells.num], [1 -1]*sum(pv)/year(), 'sat', [1 0; 0 0]);
bc    = [];
bc    = pside(bc, G, 'LEFT',  2, 'sat', [1,0], 'range', 1); 
bc    = pside(bc, G, 'RIGHT', 1, 'sat', [1,0], 'range', G.cartDims(1)); 
bf = find(any(G.faces.neighbors==0, 2));
%bc = addBC([], bf, 'pressure', G.faces.centroids(bf,1));

fluid = initSimpleFluid('mu' , [   1,  10].*centi*poise     , ...
                        'rho', [1014, 859].*kilogram/meter^3, ...
                        'n'  , [   2,   2]);
state = initResSol(G, 0, 0.0);
S     = computeMimeticIP(G, rock, 'InnerProduct', 'ip_simple');
tic
state = solveIncompFlow(state, G, S,  fluid, 'src', src, 'bc', bc, ...
                        'MatrixOutput', true);
toc
figure(1)
subplot(3,2,1:2)
   cla, for i = 1 : numel(G), plotGrid(G(i)); end, view(3)
subplot(3,2,3:4)
   plotCellData(G, state.pressure(1 : G.cells.num)); view(3)
if isfield(state, 'facePressure'),
   subplot(3,2,5)
      plot(state.facePressure);
end
subplot(3,2,6)
   plot(state.pressure);

%%

if isfield(state, 'facePressure'),
   figure(3), clf, plot(G.faces.centroids(:,1), state.facePressure, '.')
   figure(4),clf
   %fc=find(xe.facePressure >2 | xe.facePressure<1);
   [damax, fc] = max(max(state.facePressure - 2, 1 - state.facePressure));
   fc = fc(1);
   mcf = find(G.cells.faces(:,1) == fc);
end
cells = G.faces.neighbors(fc,:)
cells = cells(find(cells(:) ~= 0));
cells = unique(cells);
tcf   = find(G.cells.faces(:,1) == fc(1));

cf = mcolon(G.cells.facePos(cells), G.cells.facePos(cells+1) - 1)
plotGrid(G, cells, 'EdgeAlpha', 0.1, 'FaceColor', 'none')
plotFaces(G, fc), view(3)
