function [G]=addhybrid(G,hybrids,apertures,varargin)
% Add hybrid cells to the grid structure.
%
% A hybrid cell can be considered a lower-dimensional object in the
% geometrical grid, that still has a volume in the computational grid.
% Hybrid cells are created by converting a face into a thin cell.
% The hybrid cell consist of the same number of points as the correspoing face,
% but a volume = aperture * A assosiated with each of the hybrid cell,
% where A is the hybrid face area.
%
% In addition to the original faces (the ones that the hybrid cell was
% created from) the hybrid cells consists of hybrid faces. They are lines
% in 3d and points in 2d, but with an area = aperture * L (3d) and area =
% aperture (2d), where L is the length of the correspoing face edge. 
%
% The hybrid cells will be connected both to 'normal' cells (the ones that
% used to share a face before the the hybrid cells was squeezed in) and
% possibly other hybrid cells (if the face formed a segment of a longer
% fracture, say). The former connections are stored in G.faces.neighbors,
% for details on the ordering of faces etc see comments witin the
% subfunction 'createHybridCells' below. Connections between hybrid
% neighbors are stored in a special field G.hybridNeighbors. For flow
% calculations on the resulting grid hybrid transmissibilities between the
% hybrid faces can be computed from computeHybridTrans.m or
% computeHybridMPTrans.m
%
% If fractures (e.g. lines / planes of multiple hybrid cells) intersects, a
% small cell may be created in the intersection. However, this has a really
% small volume, and it may be desirable to eliminate it in flow
% calculations. This is done by default, unless the optional field
% addHybrid2Cells is set to true. The elimination results in non-neighbor
% cell-to-cell connections between hybrid cells that are stored in
% G.hybridNeighbors.
%
% SYNOPSIS
%
%   G = addhybrid(G,faces,apertures)
%
% PARAMETERS:
%   G         - Grid data structure.
%
%   faces     - Logical Vector of face indices. A hybrid cell is created
%               for each true face.
%
%   apertures - Vector of size G.faces.num giving the aperture of the
%               hybrid cell.
%
%
% OPTIONS:
%  addHybrid2Cells - adds hybrids cells of second kind (explanation above)
%                    and remove the internal boundary between the rest of
%                    the hybrid faces i.e no cell2cell connection are
%                    needed.  
%  addCorners      - hybrid cells of second kind are created
%                    in the corners
% EX:
%
% Create a standard grid
%       G = tensorGrid(linspace(0,10,11),linspace(0,10,11));
%       G = computeGeometry(G);
%
% Create a fracture at y=5
%       hit=find(G.faces.centroids(:,2)==5);
%
% Mark the face as a fracture face
%       G.faces.tags=zeros(G.faces.num,1);
%       G.faces.tags(hit)=1;
%
% Assign aperture
%       apt = zeros(G.faces.num,1);
%       apt(hit) = 0.001;
%
% Add the hybrid cells
%       G = addhybrid(G,G.faces.tags > 0,apt);
%
% Plot the grid
%       plotGrid(G);
%       plotFractures(G);
%
% Copyright 2011-2012 University of Bergen
%
% This file is licensed under the GNU General Public License v3.0.

opt = struct('addHybrid2Cells',false, ...
    'addCorners', false);

opt = merge_options(opt, varargin{:});

% If no hybrid faces are specified, do nothing
if ~any(hybrids)
    return
end

% apertures must be spesified to each hybrid face
assert(all(apertures(hybrids)~=0))

% creates the fracture geometry and topology and adds it to
% the grid structure.
[G,eta,hybridFaces] = createHybridCells(G,hybrids,apertures);

% the connection list between all hybrid faces are stored
% in G.hybridNeighbors.
G.hybridNeighbors = getHybridNeighbors(eta,hybridFaces);

if opt.addHybrid2Cells
    G = addHybrid2Cells(G,opt);
end

function [G,eta,hybridFaces] = createHybridCells(G,hybrids,apertures)
% creates the hybrid cells and adds them to the grid structure
%
% PARMAMETERS:
%   eta         - number of hybrid faces meeting at an intersection
%   hybridFaces - index list of all connected hybrid faces.

% number of hybrid cells
nhybrids = sum(hybrids);

% number of faces
nfaces = G.faces.num;

% maps cellfaces to cell number
cellNo = rldecode(1:G.cells.num, diff(G.cells.facePos), 2)';

% if 1 the face is hybrid, if 0 it is a normal face
G.faces.hybrid = zeros(nfaces,1);

% the cell tags stores the originating face index.
if ~isfield(G.cells,'tags')
    G.cells.tags = zeros(G.cells.num,1);
end
% if 1 the cell is hybrid, if 0 it is a normal cell
G.cells.hybrid = zeros(G.cells.num,1);

% Indices of hybrid cells, they are put in the end of the cell list
hybridCells = (G.cells.num+1:G.cells.num+nhybrids)';

% pick out the nodes of the given faces
offset = diff(G.faces.nodePos);
n = G.faces.nodes(rldecode(hybrids,offset)>0);

%hc=rldecode(find(hybrids),offset(hybrids),1);
%splittnodes(G,n,hc)

%% Copy faces

% The geometric information is kept (the hybrid elements are lower
% dimentional, they have no extension)

G.faces.nodes = [G.faces.nodes ; n];
G.faces.nodePos = [G.faces.nodePos ; G.faces.nodePos(end) + ...
    int32(cumsum(double(offset(hybrids))))];
G.faces.areas = [G.faces.areas ; G.faces.areas(hybrids)];
G.faces.normals = [G.faces.normals ; G.faces.normals(hybrids,:)];
G.faces.centroids = [G.faces.centroids ; G.faces.centroids(hybrids,:)];
G.faces.hybrid = [G.faces.hybrid ; G.faces.hybrid(hybrids)];
G.faces.tags = [G.faces.tags ; G.faces.tags(hybrids)];

% Insert hybrid cells into the neighbors in faces. The connection between
% neighbors(hybrid,1) and the hybrid cell gets the index of the old
% cell-to-cell connection, while hybrid-neighbors(hybrid,2) is stored
% towards the end of the faces, and thus have index higher than nfaces.

% Old neighbor relations, a hybrid element is about to be introduced
% inbetween
oldNeighbors = G.faces.neighbors(hybrids,:);

% Tell old cells about hybrids
G.faces.neighbors(hybrids,:) = [oldNeighbors(:,1) , hybridCells];

% And set the new hybrid relations
G.faces.neighbors = [G.faces.neighbors ; [hybridCells , oldNeighbors(:,2)]];

% Indices of hybrid faces in the original grid
hybFace=find(hybrids);

%% Then inform the cells about their new faces

% This only applies to oldNeighbors(:,2), since no update was done to the
% face in the first colmn; we simply changed the cell on the other side of
% the face

% Only the intern face index must be updated
isIntern = oldNeighbors(:,2)~=0;

% Find cells with updated faces
ind = ismember([G.cells.faces(:,1),cellNo],[hybFace(isIntern),oldNeighbors(isIntern,2)],'rows');

% Face number
fn = G.faces.num + 1 : G.faces.num + nhybrids;

% Pick out the intern faces
fn = fn(isIntern)';

% Use the same sorting for the new and old indexes
[~,I] = sort(oldNeighbors(isIntern,2));
G.cells.faces(ind,1) = fn(I);

% Update the total number of faces
G.faces.num = G.faces.num + nhybrids;

%% Connect hybrid faces

% find the mapping between fracture segments and a vertex (f)
dim = size(G.nodes.coords,2);
if(dim==2)
    faceNodes = n;
    uniqueNodes = unique(faceNodes,'rows');
    f = bsxfun(@eq,faceNodes(:,1),uniqueNodes(:,1)');
else
    i = (1:length(n))'+1;
    i(cumsum(double(offset(hybrids)))) = i(cumsum(double(offset(hybrids))))-(double(offset(hybrids)));
    faceNodes = [n,n(i)];
    faceNodes = sort(faceNodes,2);
    uniqueNodes = unique(faceNodes,'rows');
    f = bsxfun(@eq,faceNodes(:,1),uniqueNodes(:,1)') & bsxfun(@eq,faceNodes(:,2),uniqueNodes(:,2)');
end
% eta is the number of fracture segments meeting at a vertex
eta = sum(f)';

% mappings between hybrid faces and hybrid cells index
hybridFace2cell = rldecode(hybridCells,offset(hybrids));

% mapping between hybrid faces and to hybrid face index
hybridFace2face = rldecode(find(hybrids),offset(hybrids));

clear i;

% index of the vertex (i) and the faces (j).
[i,j] = find(f);

% check if we are at the boundary
tags = isBoundary(G,uniqueNodes);
isb = any(tags < 0,2);

%% fixing normals
% In the case where fractures ends within domain, the normals must be tilted
% to avoid degenerated interaction zone. 
G = fixingNormals(G,eta,isb,i,j,faceNodes,nhybrids,hybrids,offset,nfaces,hybridFace2face);


%%

% hybrid faces are created when eta>1 or when a fracture is ending
% at the boundary.
isHybridface = (isb & eta==1 | eta>1);

% maps from vertex to face
isHybridface(i) = isHybridface(j);

% pick the wanted nodes
nodes = faceNodes(isHybridface,:);

% number of hybrid faces
nhybridfaces = length(nodes);

% we store an indexlist of all the hybrid faces
% that are connected i.e. eta>1.
cs = cumsum(~isHybridface);
hybridFaces = G.faces.num+i(eta(j)>1)-cs(i(eta(j)>1));

%% Create hybrid elements

% Make hybrid faces
G.faces.nodes = [G.faces.nodes;reshape(nodes',[],1)];
G.faces.nodePos = [G.faces.nodePos;G.faces.nodePos(end)+int32(cumsum(ones(nhybridfaces,1)*(dim-1)))];
G.faces.hybrid = [G.faces.hybrid;ones(nhybridfaces,1)];
G.faces.neighbors = [G.faces.neighbors;[hybridFace2cell(isHybridface),zeros(nhybridfaces,1)]];
G.faces.centroids = [G.faces.centroids;getCentroids(G,nodes)];
a = rldecode(apertures(hybrids),offset(hybrids));
G.faces.areas = [G.faces.areas;getAreas(G,nodes,a(isHybridface))];
G.faces.normals = [G.faces.normals;getNormals(G,nodes,hybridFace2face(isHybridface))];
G.faces.tags = [G.faces.tags;zeros(nhybridfaces,1)];
G.faces.num = G.faces.num+nhybridfaces;

% Make hybrid cells

% Look through all face neighbors and pick out the ones linked to hybrid
% cells.
I1 = find(ismember(G.faces.neighbors(:,1),hybridCells));
I2 = find(ismember(G.faces.neighbors(:,2),hybridCells));
faces = [I1;I2];

% The corresponding cells
cells = [G.faces.neighbors(I1,1);G.faces.neighbors(I2,2)];

% Make sure the correct faces are linked to the cells
[cells,map] = sort(cells);
faces = faces(map);

% Update the grid structure
G.cells.faces = [G.cells.faces(:,1) ; faces]; % we dont know about the order of the faces maybe we should?
G.cells.facePos = [G.cells.facePos ; G.cells.facePos(end) + int32(find(diff(cells)));...
    G.cells.facePos(end) + length(cells)];
G.cells.centroids = [G.cells.centroids;G.faces.centroids(hybrids,:)];
G.cells.volumes = [G.cells.volumes;G.faces.areas(hybrids).*apertures(hybrids)];
G.cells.tags = [G.cells.tags;find(hybrids)];
G.cells.hybrid = [G.cells.hybrid;ones(nhybrids,1)];
G.cells.num = G.cells.num+nhybrids;

function hybridNeighbors=getHybridNeighbors(eta,hybridFaces)
% Connects the hybrid faces.
%
% PARAMETERS:
%   eta             - Number of fracture segments meeting at an vertex
%   hybridFaces     - index list of all connected hybrid faces.
%   hybridNeigbors  - Structure storing connection data for the hybrid
%                     elements. Is needed to calculated the hybrid transmissibilities.

% To or more fracture segments must meet to form a connection
eta = rldecode(eta(eta>1),eta(eta>1));

% Create the structure
hybridNeighbors.faces = [];
hybridNeighbors.n = [];
facePos(1) = 1;
neighbors = [];

% Find the connections.
startIndex = 0;
for k = unique(eta)'
    
    % pick out the hybrid faces for a given eta
    faces = hybridFaces(eta==k);
    hybridNeighbors.faces = [hybridNeighbors.faces;faces];
    
    % there will be sum_i^(eta-1) i connections for each vertex
    n1 = rldecode((1:k)',repmat(k-1:-1:0,1));
    n2 = n1+mcolon(ones(1,k-1),k-1:-1:1)';
    nc = length(faces)/k;
    n = sum(1:k-1);
    index = rldecode((0:k:length(faces)-1)',n) + startIndex;
    neighbors = [neighbors; [repmat(n1,nc,1) + index,repmat(n2,nc,1) + index]]; %#ok<AGROW>
    startIndex = neighbors(end,2);
    
    % store the index such that faces(facePos(i):facePos(i+1)-1) return all
    % faces meeting at a vertex i
    facePos = [facePos ; facePos(end) + (k : k : length(faces))']; %#ok<AGROW>
    
    % Number of connections.
    hybridNeighbors.n = [hybridNeighbors.n ; n * ones(nc,1)];
    
end

% update the structure.
hybridNeighbors.neighbors = neighbors;
hybridNeighbors.facePos = facePos;

function G = addHybrid2Cells(G,opt)
% Adds hybrid cells of secound kind and remove interanl
% boundaries between the hybrid faces.

% numbers of hybrid segments meeting at a vertex
eta = diff(G.hybridNeighbors.facePos);
eta2 = rldecode(eta,eta);

% the hybrid nodes
nodes = G.faces.nodes(G.faces.nodePos(G.hybridNeighbors.faces));

% hybrid 2 cells are created only when more then two
% hybrid sements meet at a vertex or if wanted at the corners
if opt.addCorners
    isCorner = returnCorners(G,nodes);
    is_hybfaces = eta2 > 2 | isCorner;
else
    is_hybfaces = eta2 > 2;
end

% faces not neighbouring hybrid cells
faces = G.hybridNeighbors.faces(~is_hybfaces);

% faces neighbouring hybrid cells
hyb_faces = G.hybridNeighbors.faces(is_hybfaces);

% number of neighboring faces
hyb_eta = eta(is_hybfaces(G.hybridNeighbors.facePos(1:end-1)));

% number of hybrid cells
num_hc = numel(hyb_eta);

% connect faces to the hybrid cells
j = rldecode(1:num_hc,hyb_eta,2)';
G.faces.neighbors(hyb_faces,2) = j+G.cells.num;

% create the hybrid cells.
G.cells.facePos = [G.cells.facePos; int32(cumsum(hyb_eta))+G.cells.facePos(end)];
G.cells.faces = [G.cells.faces; hyb_faces];
G.cells.num = G.cells.num + num_hc;
G.cells.centroids = [G.cells.centroids; G.faces.centroids([hyb_faces(find(diff(j))) ...
                                                         ; hyb_faces(end)],:)];
G.cells.tags = [G.cells.tags; zeros(num_hc,1)];
G.cells.hybrid = [G.cells.hybrid; 2*ones(num_hc,1)];
G.cells.volumes = [G.cells.volumes; zeros(num_hc,1)];

% remove the internal boundary between the rest of the
% hybrid faces.
N = reshape(faces,2,[])';
G = removeInternalBoundary(G,N);

%% helpers
function G = fixingNormals(G,eta,isb,i,j,faceNodes,nhybrids,hybrids,offset,nfaces,hybridFace2face)
% When a fracture ends, special treatment is needed of the normal vectors
% to avoid degenerated interaction regions. Do this by perturbing the
% normal vectors slightly.
% An alternative approach, which is more consistent with an equidimensional
% representation is to use normal vectors associated with the expanded
% fracture. With the uncertainty in fracture descriptions in mind, we have
% not found this worthwhile.

isHybridfaceEta1 = (~isb & eta==1);
hybridFace2thereFace = rldecode((1:nhybrids)',offset(hybrids))+nfaces-nhybrids;
isHybridfaceEta1(i) = isHybridfaceEta1(j);
nn = -G.faces.centroids(hybridFace2face(isHybridfaceEta1),:)+getCentroids(G,faceNodes(isHybridfaceEta1,:));

% the normals is tilted eps towards the vertex.
nn = bsxfun(@times,nn,sqrt(eps)./sqrt(sum(nn.^2,2)));

nHere  = G.faces.normals(hybridFace2face(isHybridfaceEta1),:);
nThere = G.faces.normals(hybridFace2thereFace(isHybridfaceEta1),:);
nHere  = bsxfun(@times,nHere,1./sqrt(sum(nHere.^2,2)));
nThere = bsxfun(@times,nThere,1./sqrt(sum(nThere.^2,2)));

nHere = nHere+nn;
nThere = nThere+nn;

nHere = bsxfun(@times,nHere,G.faces.areas(hybridFace2face(isHybridfaceEta1))./sqrt(sum(nHere.^2,2)));
nThere = bsxfun(@times,nThere,G.faces.areas(hybridFace2thereFace(isHybridfaceEta1))./sqrt(sum(nThere.^2,2)));
G.faces.normals(hybridFace2face(isHybridfaceEta1),:) = nHere;
G.faces.normals(hybridFace2thereFace(isHybridfaceEta1),:) = nThere;


function c = getCentroids(G,nodes)
% return the centroid of the given nodes

dim = size(G.nodes.coords,2);
if(dim==2)
    c = G.nodes.coords(nodes,:);
    return
end
v = G.nodes.coords(nodes(:,1),:)-G.nodes.coords(nodes(:,2),:);
c = G.nodes.coords(nodes(:,2),:)+v./2;

function n = getNormals(G,nodes,faces)
% returns the normal given nodes defining the line/plane
dim = size(G.nodes.coords,2);

% for 2d it is straigh forward
if(dim==2)
    n = G.nodes.coords(nodes,:)-G.faces.centroids(faces,:);
    areas = G.faces.areas(G.faces.hybrid==1);
    n = bsxfun(@times,n,areas./sqrt(sum(n.*n,2)));
    return
end

% the normal is given by the cross product between the face normal and
% a vector going between the end points of the hybrid face

% a vector going between the end points of the hybrid face
v = G.nodes.coords(nodes(:,2),:)-G.nodes.coords(nodes(:,1),:);

% If the hybrid face has a neighbor, get its normal
neighbors = G.faces.neighbors(G.faces.hybrid==1,:);
i = neighbors(:,1)~=0;
j = neighbors(:,2)~=0;
np(i,:) = G.faces.normals(faces,:);

% Cross them to get the normal
n = cross(np,v);

% Then we have to make sure it points in the right direction.

% we define a vector between the center of the adjesant faces or
% between the center of the hybrid face and one adjesant face in the
% case where the hybrid face is at the bondary.

x1 = G.faces.centroids(G.faces.hybrid==1,:);
x2 = x1;
x1(i,:) = G.faces.centroids(faces(i),:);
x2(j,:) = G.faces.centroids(faces(j),:);
r = x2-x1;

% we want the normal pointing outward
signchange = dot(r,n,2)<0;
n(signchange,:) = -n(signchange,:);

% Then we make sure that the length of the normals equals the area.
areas = G.faces.areas(G.faces.hybrid==1);
n = bsxfun(@times,n,areas./sqrt(sum(n.*n,2)));

function A = getAreas(G,n,a)
% returns the hybrid area, given some nodes (n) and the aperture (a).

dim = size(G.nodes.coords,2);
if(dim==2)
    A = a;
    return
end
d = G.nodes.coords(n(:,1),:)-G.nodes.coords(n(:,2),:);
A = sqrt(sum(d.^2,2)).*a;

function tags = isBoundary(G,nodes)
% return boundary tags (-1 -- -6 ) if the nodes lays on the boundary.
% For the current implementation, the only information needed is whether
% the tag is less than zero or not.

minc = min(G.nodes.coords);
maxc = max(G.nodes.coords);
dim = size(G.nodes.coords,2);
tags = zeros(size(nodes,1),dim*2);

i = all(reshape(G.nodes.coords(nodes,1)==minc(1),size(nodes)),2);
tags(i,1) = -1;

i = all(reshape(G.nodes.coords(nodes,1)==maxc(1),size(nodes)),2);
tags(i,2) = -2;

i = all(reshape(G.nodes.coords(nodes,2)==minc(2),size(nodes)),2);
tags(i,3) = -3;

i = all(reshape(G.nodes.coords(nodes,2)==maxc(2),size(nodes)),2);
tags(i,4) = -4;

if dim == 2
    return
end

i = all(reshape(G.nodes.coords(nodes,3)==minc(3),size(nodes)),2);
tags(i,5) = -5;

i = all(reshape(G.nodes.coords(nodes,3)==maxc(3),size(nodes)),2);
tags(i,6) = -6;



function isCorner = returnCorners(G,nodes)
box = [min(G.nodes.coords),max(G.nodes.coords)]; %x1 y1 x2 y2 

c1 = isEq(G.nodes.coords(nodes,1),box(1)) &  isEq(G.nodes.coords(nodes,2),box(2));
c2 = isEq(G.nodes.coords(nodes,1),box(1)) &  isEq(G.nodes.coords(nodes,2),box(4));
c3 = isEq(G.nodes.coords(nodes,1),box(3)) &  isEq(G.nodes.coords(nodes,2),box(2));
c4 = isEq(G.nodes.coords(nodes,1),box(3)) &  isEq(G.nodes.coords(nodes,2),box(4));
isCorner = c1 | c2 | c3 | c4;

function i = isEq(a,b)
% To numbers (coordinates) are defined as equal if they differ with less
% than sqrt(eps).
i = abs(a-b) < sqrt(eps);
