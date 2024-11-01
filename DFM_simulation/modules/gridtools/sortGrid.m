function G = sortGrid(G)
%Permute nodes, faces and cells to sorted form
%
% SYNOPSIS:
%   G = sortGrid(G)
%
% PARAMETERS:
%   G       - Grid structure as described by grid_structure.
%
% RETURNS:
%   G       - Modified grid structure with nodes in natural order, sorted
%             z-, y- and x-components, faces sorted by corner node numbers
%             and cells sorted by face numbers.
%
% COMMENTS:
%
%  Note that the ordering of faces and cells is NOT the natural ordering
%  as generated by cartGrid and tensorGrid.   
%
% SEE ALSO:
%  compareGrids

%{
Copyright 2009, 2010, 2011, 2012 SINTEF ICT, Applied Mathematics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

% $Date: 2012-12-19 16:47:24 +0100 (Wed, 19 Dec 2012) $
% $Revision: 10500 $

   if ~isfield(G.cells, 'facePos'),
       G.cells.facesPos = int32(cumsum([1; double(G.cells.numFaces)]));
       G.cells = rmfield(G.cells, 'numFaces');
   end
   
   if ~isfield(G.faces, 'nodePos'),
       G.faces.nodePos = int32(cumsum([1; double(G.faces.numNodes)]));
       G.faces = rmfield(G.faces, 'numNodes');
   end
   
   G = sortNodes(G);
   G = sortFaces(G);
   G = sortCells(G);
   G.type = [G.type, mfilename];
end  

function g = sortNodes(g)
   n              = fliplr(g.nodes.coords);
   [n, map]       = sortrows(n);
   g.nodes.coords = fliplr(n);
   inversemap     = accumarray(map, 1:g.nodes.num);
   g.faces.nodes  = inversemap(g.faces.nodes);
end

function g = sortFaces(g)
   
   % Start each faceNodes with smallest node index
   e            = expand(g.faces.nodePos, g.faces.nodes);

   if any(diff(sort(e))==0),
      error('There are duplicate nodes in one or more faces');
   end
   [r, dummy]   = find(bsxfun(@eq, min(e, [], 1), e));%#ok
   ix           = rot(g.faces.nodePos, r-1);
   g.faces.nodes  = g.faces.nodes(ix);
   e            = expand(g.faces.nodePos, g.faces.nodes);
  
   % sort faces
   
   [E, I]=  sortrows(e');
   
   E = E';
   invI = accumarray(I, 1:g.faces.num);
   g.faces.nodes = E(~isnan(E));
   g.faces.neighbors = g.faces.neighbors(I, :);
   numNodes = diff(g.faces.nodePos);
   g.faces.nodePos = int32(cumsum([1; double(numNodes (I))]));
   g.cells.faces(:,1) = invI(g.cells.faces(:,1));
 
end

function g = sortCells(g)
   % Start each cellFaces with smallest node index
   e            = expand(g.cells.facePos, g.cells.faces(:,1));
   [r, dummy]   = find(bsxfun(@eq, min(e(:,:,1), [], 1), e(:,:,1)));%#ok
   if numel(r) ~= size(e,2),
      error('There are duplicate faces in one or more cells');
   end
   ix           = rot(g.cells.facePos, r-1);
   g.cells.faces(:,1)  = g.cells.faces(ix,1);
   
   % Primitive processing may skip face tags
   if size(g.cells.faces, 2) == 2,
      g.cells.faces(:,2)  = g.cells.faces(ix,2);
   end
   e            = expand(g.cells.facePos, g.cells.faces(:,1));
   e = sort(e);
   [E, I]=  sortrows(e');
   E = E';
   invI = accumarray(I, 1:g.cells.num);
   g.cells.faces(:,1) = E(~isnan(E));
   i = g.faces.neighbors ~= 0;
   g.faces.neighbors(i) = invI(g.faces.neighbors(i));
   numFaces = diff(g.cells.facePos);
   numFaces = numFaces(I);
   g.cells.facePos = int32(cumsum(double([1;numFaces])));
   if isfield(g.cells, 'indexMap'),
      g.cells.indexMap = g.cells.indexMap(I);
   end
   % cellFaces(:,2) buggy
end


function ix = rot(pos, r)
   v  = expand(pos, 1:pos(end)-1);
   j  = mod(bsxfun(@plus,  (0:size(v, 1)-1)', r(:)'), size(v, 1))+1;
   i  = repmat((0:size(v, 2)-1)*size(v, 1), [size(v, 1), 1]);
   ix = compress(v(i+j));
end


function [v, ix] = expand(pos, value)
   num = double(diff(pos));
   
   m   = numel(num);
   n   = max(num);
   off = reshape((0 : m - 1) .* n, [], 1);
   ix  = mcolon(off + 1, off + num);
   v     = nan([n, m]);
   v(ix) = value;
end
   
function [value, pos] = compress(v)
   num   = sum(~isnan(v));
   pos   = cumsum([1, num])';
   value = v(~isnan(v));
end
   


