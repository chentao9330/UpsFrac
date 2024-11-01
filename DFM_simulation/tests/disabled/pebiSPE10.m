function H=pebiSPE10
addpath([ROOTDIR,'projects/spe10/'])
make_spe10_data
   figure(2),clf,
   H = makeGrid(1);
   H = computeGeometry(H);H.cells.num
   
   rock     = SPE10_rock(1);
   rock.perm = rock.perm(:,1:2)*milli*darcy;
   rock.poro = max(0.001,rock.poro);
   i         = max(1, round(H.cells.centroids(:,1)/(20*ft)));
   j         = max(1, round(H.cells.centroids(:,2)/(10*ft)));
   rock.perm = rock.perm(sub2ind([60, 220], i, j));
   rock.poro = rock.poro(sub2ind([60, 220], i, j));
   
if false
   wellpos = {[312.2863   12.7446], ...
              [269.1509  640.1691]};
   %src = [];
   %src = addSource(src, findEnclosingCell(H, wellpos{1}),  1e-2, 'sat', 1);
   %src = addSource(src, findEnclosingCell(H, wellpos{2}), -1e-2, 'sat', 0);
   
   W =[];
   W = addWell(W, H, rock, findEnclosingCell(H, wellpos{1}), 'Type','bhp', 'val',500*barsa);
   W = addWell(W, H, rock, findEnclosingCell(H, wellpos{2}), 'Type','bhp', 'val',100*barsa);
else
   % Use wells (productivity index, position, pressure) from the original 
   % case (simulate_spe10)

   g          = computeGeometry(cartGrid([60,220], [60*20*ft, 220*10*ft]));
   cells      = sub2ind([60, 220], ...
                        [1, 60,  60,   1,  30], ...
                        [1,  1, 220, 220, 110])';
   pos        = g.cells.centroids(cells,:);
   wellpos    = num2cell(pos, 2);
   
   WI = [...
      0.034449455252505, ...
      0.270623680030513, ...
      0.013093819440371, ...
      0.002821113544878, ...
      0.001759001155289]*darcy;
   W =[];
   W = addWell(W, H, rock, findEnclosingCell(H, wellpos{1}), 'Type','bhp', 'val', 4000*psia, 'WI', WI(1));
   W = addWell(W, H, rock, findEnclosingCell(H, wellpos{2}), 'Type','bhp', 'val', 4000*psia, 'WI', WI(2));
   W = addWell(W, H, rock, findEnclosingCell(H, wellpos{3}), 'Type','bhp', 'val', 4000*psia, 'WI', WI(3));
   W = addWell(W, H, rock, findEnclosingCell(H, wellpos{4}), 'Type','bhp', 'val', 4000*psia, 'WI', WI(4));
   W = addWell(W, H, rock, findEnclosingCell(H, wellpos{5}), 'Type','bhp', 'val',10000*psia, 'WI', WI(5));
   
   clear g pos

end
   fluid = initSimpleFluid('mu', [0.3,3]*centi*poise, 'rho', [1000, 10],'n',[2,2]);
   
   gravity on
   gravity([0,-10]);
   gravity off;
   
   IP = computeMimeticIP(H, rock);
   x  = initResSol(H, 0, 0);
   
   
   subplot(1, 2, 1);
   plotCellData(H, log10(rock.perm/darcy));colorbar,axis equal tight off
   for i=1:numel(wellpos),
      hold on;
      plot(wellpos{i}(1), wellpos{i}(2), 'o',...
           'markerfacecolor', 'w', 'markersize',5);
   end
   box on
   drawnow;
   
   subplot(1,2,2);
   cla;
   for i=1:numel(wellpos),
      hold on;
      plot(wellpos{i}(1), wellpos{i}(2), 'o',...
           'markerfacecolor', 'w', 'markersize',5);
   end
   
   for i=1:50,
      x  = solveIncompFlow(x, H, IP, fluid, 'wells',W);
      x = implicitTransport(x, H, 100*day, rock, fluid, 'wells', W,'verbose',true);
      plotCellData(H, x.s(:,1));axis equal tight off
      box on
      drawnow;
   end

end


function c = findEnclosingCell(G, pt)
%
% Find cells with nearest centroid... Then check carefully.
%
   
   %G     = sortEdges(G);
   edges = reshape(G.faces.nodes, 2, [])';
   
   % For each edge, check which side x is.
   x = G.nodes.coords(:,1); 
   y = G.nodes.coords(:,2);
   
   a = [-diff(y(edges), 1, 2), diff(x(edges), 1, 2)];
   b = bsxfun(@minus, pt, [x(edges(:,1)), y(edges(:,1))]);
   
   v      = sum(a.*b, 2);
   cellno = rldecode(1:G.cells.num, diff(G.cells.facePos), 2)';
   sgn    = -1 + 2*(cellno==G.faces.neighbors(G.cells.faces(:,1), 1));
   V      = sgn.*v(G.cells.faces(:,1))>=0;
   
   i = accumarray(cellno, V, [G.cells.num, 1], @(x) all(x));
   c = find(i);
end


function G =  makeGrid(type)
   %points = [0,0; 60,0; 60,220; 0,220; 10,50; 50,100];%meter
   points = [0,0; 60,0; 60,220; 0,220];%meter
   points = bsxfun(@times, points, [20*ft, 10*ft]);
   edges  = [1,2;2,3;3,4;4,1]
   G = triangleGrid(points, edges, 'maxArea', 30);
   G.nodes.coords = G.nodes.coords;
   G = pebi(G)
end

