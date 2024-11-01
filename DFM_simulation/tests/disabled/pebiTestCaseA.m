function pebiTestCaseA
   %H =  cartGrid([50,200],[60,220]);%
   H = makeGrid(1);
   H = computeGeometry(H);
   
   wells = {[52.5,8.1], [31.5,215.2]};
   wells = {[182.8800, 42.1551], [171.1158, 634.2870]};
   src = [];
   src = addSource(src, findEnclosingCell(H, wells{1}),  1e-2, 'sat', 1);
   src = addSource(src, findEnclosingCell(H, wells{2}), -1e-2, 'sat', 0);
  
   load Udata;
   kx        = squeeze(KU(1,:,:,1))*0+10*milli*darcy;
   phi       = squeeze(pU(:,:,1))*0+0.3;
   i         = max(1, round(H.cells.centroids(:,1)/(20*ft)));
   j         = max(1, round(H.cells.centroids(:,2)/(10*ft)));
   rock.perm = kx(sub2ind(size(kx), i, j));
   %rock.perm = repmat(1*micro*darcy, [H.cells.num, 1]);
   rock.poro = max(1e-2, phi(sub2ind(size(kx), i, j)));
   
   fluid = initSimpleFluid('mu', [1,1]*centi*poise, 'rho', [1000, 10],'n',[2,2]);
   
   gravity on
   gravity([0,-10]);
   gravity off;
   
   IP = computeMimeticIP(H, rock);
   x  = initResSol(H, 0, 0);
   
   clf;
   subplot(1, 2, 1);cla
   plotCellData(H, log10(rock.perm));axis equal tight off
   for i=1:numel(wells),
      hold on;
      plot(wells{i}(1), wells{i}(2), 'o',...
           'markerfacecolor', 'w', 'markersize',10);
   end
   box on
   shading faceted
   a = axis; 
   e = 0.01*(a(2)-a(1));
   f = 0.01*(a(4)-a(3));
   a = [a(1)-e, a(2)+e, a(3)-f,a(4)+f];
   axis(a);
   drawnow;
   
   subplot(1,2,2);
   cla;
   
   for i=1:57,
      x  = solveIncompFlow(x, H, IP, fluid, 'src', src);%, ...
      %     'LinSolve',@(A,b) agmg(A,b));
      x = explicitTransportNew(x, H, 1*day, rock, fluid, 'src',src,'verbose',false);
      plotCellData(H, x.s(:,1));axis equal tight off
      title(sprintf('%d days', i));
      box on
      axis(a);
      drawnow;
   end
%print -dpng -r300 /home/jrn/devel/c++/trunk /reorder/testcases/2d-pebi/homogeneous-2079/pebi-2079.png
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
   points = [0,0; 60,0; 60,220; 0,220; 10,50; 50,100];%meter
   points = bsxfun(@times, points, [20*ft, 10*ft]);
   edges  = [1,2;2,3;3,4;4,1;5,6]
   G = triangleGrid(points, edges, 'maxArea', 1000);
   G.nodes.coords = G.nodes.coords;
   G = pebi(G)
end

