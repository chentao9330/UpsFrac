function gravitySplitting()
   % Example relies on 'streamlines' module.  Capture module state to avoid
   % spamming controlling enviroment.
   MODULES = mrstModule;
   mrstModule add streamlines

   G = cartGrid([50,50]);

   G = computeGeometry(G);
   x = initResSol(G, 0);
   x.s(G.cells.centroids(:,2)>25) = 0.5;
   
   rotateNodes = @(p, t) p*[cos(t),sin(t);-sin(t), cos(t)];
   G.nodes.coords = rotateNodes(G.nodes.coords, 210*pi/180);
   
   G = computeGeometry(G);
   
   
   try
      rock = SPE10_rock(1:50,1:50,40);
      rock.perm = convertFrom(rock.perm(:,1), milli*darcy);
      rock.poro = max(1e-4, rock.poro);
   catch ME %#ok
      rock.perm = ones(G.cells.num, 1)*100*milli*darcy;
      rock.poro = ones(G.cells.num, 1);
   end
   
   fluid = initSimpleFluid('mu',  [1,      1]*centi*poise,...
      'rho', [1000, 500],            ...
      'n',   [2,      2]);
   
   gravity([0, 10])
   
   
   
   
   T  = computeTrans(G, rock);
   IP = computeMimeticIP(G, rock);
   
   
   plotGrid(G);
   
   strength = 1e1;
   src = [];
   src = addSource(src, 1,            strength*meter^3/day, 'sat', [1,0]);
   src = addSource(src, G.cells.num, -strength*meter^3/day, 'sat', [1,0]);
   
   
   %psolve = @(state, dt) ...
   %         incompTPFA(state, G, T, fluid,'src',src);
   psolve = @(state, dt) ...
      solveIncompFlow(state, G, IP, fluid,'src',src);
   tsolvea = @(state, dt) ...
      implicitTransport(state, G, dt, rock, fluid, ...
      'src', src, 'Verbose', true);
   tsolvec = @(state, dt) ...
      implicitTransport(state, G, dt, rock, fluid, ...
      'Verbose', true);
   clf,
   t = 0;
   for step = 1 : 50,
      dt = 100*day;
      
      gravity off;
      xa     = psolve(x, dt);
      
      gravity on
      x     = psolve(x, dt);
      
      aflux = xa.flux;
      gflux = x.flux;
      
      if true,
         
         % Advection
         disp('Advection');
         x.flux = aflux;
         [x,r] = tsolvea(x, dt);
         
         % Convection
         disp('Convection');
         gravity off
         x.flux = gflux - aflux;
         [x,r] = tsolvec(x, dt);
         
         % Segregation
         disp('Segregation');
         gravity on;
         x.flux = zeros(G.faces.num, 1);
         [x,r] = tsolvec(x, dt);
      else
         x.flux = gflux;
         [x,r] = tsolve(x, dt);
      end
      
      t = t + dt;
      fprintf('Transport step: %s\n', r.str);
      
      subplot(1,3,1)
      plotCellData(G, x.s(:,1));
      view(2);
      axis equal tight off
      
      %[ii,jj] = ndgrid(3:5:G.cartDims(1), 3:5:G.cartDims(2));
      %slcells = reshape(sub2ind(G.cartDims, ii,jj), [], 1);
      slcells = (G.cartDims(1):G.cartDims(1)-1:G.cells.num)';
      subplot(1,3,2);
      cla,
      %plotGrid(G, 'facea',0.05, 'edgea',0.05);
      x.flux = gflux-aflux;
      streamline(pollock(G, x, slcells, 'maxsteps',800));
      view(2);
      axis equal tight
      box on
      set(gca, 'xtick',[],'ytick',[]);
      
      subplot(1,3,3);
      cla,
      %plotGrid(G, 'facea',0.05, 'edgea',0.05);
      x.flux = aflux;
      streamline(pollock(G, x, slcells, 'maxsteps',800));
      streamline(pollock(G, x, slcells, 'maxsteps',800, 'reverse', true));
      view(2);
      axis equal tight; a=axis;a=a+0.05*[-1,1,-1,1];axis(a);
      box on
      set(gca, 'xtick',[],'ytick',[]);
      
      drawnow
   end

   % Restore caller's module state (effectively remove 'streamlines' module
   % if we activated the module).
   mrstModule clear
   mrstModule('add', MODULES{:})
end

function rock = SPE10_rock(varargin)
% SPE10_rock -- Define rock properties for Model 2 of tenth SPE comparative
%               solutions project.
%
% SYNOPSIS:
%   rock = SPE10_rock
%   rock = SPE10_rock(layers);
%
% PARAMETERS:
%   layers - Which of the 85 model layers to include in a specific test.
%            OPTIONAL.  If unspecified, all 85 layers (for a total of
%            60-by-220-by-85 == 1,122,000 cells) are included.
%
%            Some possible values are
%               layers = ( 1 : 35).';  %  Tarbert formation
%               layers = (36 : 85).';  %  Upper Ness formation
%
% RETURNS:
%   rock - Rock structure having fields 'perm' and 'poros' pertaining to
%          the specified layers.
%
% EXAMPLE:
%   rock = SPE10_rock(85)
%
% SEE ALSO:
%   SPE10_setup.

   %% Define layers
   layers = (1 : 85);
   if nargin > 0 && isnumeric(varargin{1}),
      iindex = 1:60;
      jindex = 1:220;
      layers = varargin{1};
   end
   if nargin > 2 && isnumeric(varargin{1}),
      iindex = varargin{1};
      jindex = varargin{2};
      layers = varargin{3};
   end
   
   
   %% Retrieve data from file
   data = load(fullfile(ROOTDIR, 'projects', 'em-2008', 'spe10_rock'));
   rock = data.rock;
   
   %% Extract only required layers from data
   layerIx = myindex([60, 220, 85], iindex, jindex, layers);
   rock.perm = rock.perm(layerIx(:), :);
   rock.poro = rock.poro(layerIx(:));
end


function k = myindex(sz, i, j, k)
      [I,J,K] = ndgrid(i, j, k);
      k       = reshape(sub2ind(sz, I, J, K), [], 1);
end
