function [G, Gt, rock, rock2D, bcIxVE] = makeJohansenVEgrid()
%Make an VE model based upon a data set of the Johansen formation
%
% SYNOPSIS:
%  [G, Gt, bcIx, bcIxVE, rock, rock2D] = makeJohansenVEgrid()
%
% PARAMETERS:
%   G      - Data structure for 3D grid
%   Gt     - Data structure for topsurface grid
%   rock   - Data structure for 3D rock parameters
%   rock2D - Data structure for rock parameters for topsurface grid
%   bcIxVE - Index for pressure boundary conditions in topsurface grid
%
% DESCRIPTION:
%  The Johansen formation is a candidate site for large-scale CO2 storage
%  offshore the south-west coast of Norway. We consider a sector model that
%  has been constructed based on available seismic and well data and stored
%  in the Eclipse input format (GRDECL). The model has five vertical layers
%  in the Johansen formation and five shale layers above and one below in
%  the Dunhil and Amundsen formations. The shale layers are removed and we
%  construct a 2D grid of the top surface, assuming that the major
%  fault is sealing. Moreover, we identify all outer boundaries that are
%  open to flow. A more thorough presentation of the geological model can be
%  found in the script <matlab:edit('showJohansen.m') showJohansen.m>
%
%  The grid and rock structured constructed in the following can be used for
%  subsequent 3D and/or VE simulations of CO2 injection and migration and
%  are therefore stored to file to avoid time-consuming processing.
%
%  The data files necessary to run the example can be downloaded from the
%  <http://www.sintef.no/Projectweb/MatMorA/Downloads/Johansen/ MatMoRA
%  website>.
%
% SEE ALSO:
%   runJohansenVE

try
   disp(' -> Reading Johansen.mat');
   d = fileparts(mfilename('fullpath'));
   load(fullfile(d,'Johansen'));
   return;
catch me
   disp(' -> Reading failed, constructing grid models');
end

%% Load model and construct VE grid
% Load grid geometry - you will most likely have to change the path,
% depending upon where you have stored the Johansen data-set
try
   sector = fullfile(ROOTDIR, 'examples', 'data', 'johansen', 'NPD5');
   fprintf('    Reading %s\n', sector);
   grdecl = readGRDECL([sector, '.grdecl']);
catch me
  disp('Download data from:')
  disp('http://www.sintef.no/Projectweb/MatMorA/Downloads/Johansen/')
  disp(['Put data in ', fullfile(ROOTDIR, 'examples', 'data', 'johansen')]);
end

% Load permeability and porosity
K = reshape(load([sector, '_Permeability.txt'])', [], 1);
p = reshape(load([sector, '_Porosity.txt'])',     [], 1);
grdecl.PERMX=K;grdecl.PERMY=K;grdecl.PERMZ=0.1*K;
grdecl.PORO=p;
clear sector;

% Remove shale layers Dunhil and Amundsen
disp(' -> Constructing grids ...');
grdecl.ACTNUM(K.*milli*darcy<0.11 *milli*darcy) = 0;

% Construct grid structure. By removing Dunhil and Amundsen, the grid will
% consist of multiple components. We choose the largest one as our
% reservoir.
G = processGRDECL(grdecl);
G = computeGeometry(G(1));
clear grdecl

% Construct VE grid
[Gt, G] = topSurfaceGrid(G);

% Construct structure with petrophyiscal data. Vertical permeability is set
% to 0.1 times the horizontal permeability. NB!
rock.perm = bsxfun(@times, [1 1 0.1], K(G.cells.indexMap)).*milli*darcy;
rock.poro = p(G.cells.indexMap);
rock2D    = averageRock(rock, Gt);
clear p K;

%% FIND FAULT - set inner boundary
% The main fault is assumed to be sealing and must therefore be represented
% as an inner boundary in the 2D grid.  Start by locating cells on each
% side of fault, which is defined to be between cells with index i=43 and
% i=44 for j<44.
cells2D_1 = find(Gt.cells.ij(:,1) == 43 & Gt.cells.ij(:,2) <= 44);
cells2D_2 = find(Gt.cells.ij(:,1) == 44 & Gt.cells.ij(:,2) <= 44);

% Plot the cells on opposite sides
figure;
plotGrid(Gt, 'faceColor', 'none');
plotGrid(Gt, cells2D_1, 'faceColor', 'r')
plotGrid(Gt, cells2D_2, 'faceColor', 'g')
axis tight off,
title('Cells on opposite sides of sealing fault'),

% Find the faces at the fault. Construct a mapping facesMat defined such
% that facesMat(i,j)=k if face <k> is shared by cells <i> and <j>.
facesMat = sparse(double(Gt.faces.neighbors(:,1))+1, ...
   double(Gt.faces.neighbors(:,2))+1, 1:Gt.faces.num);
facesMat = facesMat + facesMat';
faultFaces2D = diag( facesMat(cells2D_1+1, cells2D_2+1) );

% Make internal boundary and compute the geometry of the resulting grid
% Gt = makeInternalBoundary(Gt, faultFaces2D);
% Gt = computeGeometryVE(Gt);

% The function 'topSurfaceGrid' does not handle faults correctly when the
% cells on opposite sides are not phyiscally in contact with each other.
% Instead of producing a discontinuity in the 'z' value, an average value
% is used.  Hence, we need to manually reset the 'z' value of these cells
% (marked in red and green in the plot) to avoid an incorrect flow
% Gt.cells.z([cells2D_1; cells2D_2]) = ...
%    G.faces.centroids(Gt.cells.map3DFace([cells2D_1; cells2D_2]), 3);
clear cells2D_1 cells2D_2 facesMat faultFaces2D


%% FIND PRESSURE BOUNDARY
% Setting boundary conditions is unfortunately a manual process and may
% require some fiddling with indices, as shown in the code below. Here, we
% identify the part of the outer boundary that is open, i.e., not in
% contact with one of the shales (Dunhil or Amundsen).

% boundary 3D
nx = G.cartDims(1); ny=G.cartDims(2); nz=G.cartDims(3);
ix1 = boundaryFaceIndices(G, 'BACK', 1:nx-6, 1:4, 1:nz);
ix2 = boundaryFaceIndices(G, 'LEFT', 1:20,   1:ny, 1:nz);
ix3 = boundaryFaceIndices(G, 'RIGHT', 1:nx, ny-10:ny, 1:nz);
ix4 = boundaryFaceIndices(G, 'FRONT', 1:nx/2-8, ny/2:ny, 1:nz);

figure;
subplot(1,2,1)
plotGrid(G, 'faceColor', 'none', 'EdgeAlpha', 0.1)
plotFaces(G, ix1, 'r');
plotFaces(G, ix2, 'g');
plotFaces(G, ix3, 'y')
plotFaces(G, ix4, 'm');
view(-7,60), axis tight off;
title('3D grid');

% boundary 2D
nx = Gt.cartDims(1); ny=Gt.cartDims(2);
ix1 = boundaryFaceIndices(Gt, 'BACK',  1:nx-6, 1:4, []);
ix2 = boundaryFaceIndices(Gt, 'LEFT',  1:20, 1:ny,  []);
ix3 = boundaryFaceIndices(Gt, 'RIGHT', 1:nx, ny-10:ny, []);
ix4 = boundaryFaceIndices(Gt, 'FRONT', 1:nx/2-8, ny/2:ny, []);

%remove faces connected to main fault
ix1 = ix1(Gt.faces.centroids(ix1,2)>6.714*1e6);
ix2 = ix2(Gt.faces.centroids(ix2,1)>5.4*1e5);
%
subplot(1,2,2)
plotGrid(Gt, 'faceColor', 'none', 'EdgeAlpha', 0.1)
plotGrid(Gt, sum(Gt.faces.neighbors(ix1,:),2), 'faceColor', 'r')
plotGrid(Gt, sum(Gt.faces.neighbors(ix2,:),2), 'faceColor', 'g')
plotGrid(Gt, sum(Gt.faces.neighbors(ix3,:),2), 'faceColor', 'y')
plotGrid(Gt, sum(Gt.faces.neighbors(ix4,:),2), 'faceColor', 'm')
axis tight off
title('2D grid of top surface');

bcIxVE = [ix1; ix2; ix3; ix4];
clear ix1 ix2 ix3 ix4 nx ny nz

%% Store data
disp(' -> Writing Johansen.mat')
save(fullfile(d,'Johansen'), 'G', 'Gt', 'rock', 'rock2D', 'bcIxVE');
end
