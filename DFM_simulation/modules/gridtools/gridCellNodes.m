function [n, pos] = gridCellNodes(G, c)
% Find nodes corresponding to a set of cells
%
% SYNOPSIS:
%   [n, pos] = gridCellNodes(G, c)
%
% PARAMETERS:
%   G    - Grid structure
%   c    - Cells where the fine nodes are desired
%
%
% RETURNS:
%   pos  - indirectionmap into n. The nodes of cell c(i) is found at
%   positions p(i):p(i+1)-1 in n
%   n    - node positions in G
%

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
    
    % Find number of faces per cell
    nf = diff([G.cells.facePos(c), G.cells.facePos(c+1)], [],2);
    
    % Find the cell index of each face
    cellno = rldecode(1:numel(c), nf, 2) .';
    
    % Number of nodes per face
    nnode = diff(G.faces.nodePos);
    
    % Find faces of cell subset
    cf = G.cells.faces(mcolon(G.cells.facePos( c ), ...
                              G.cells.facePos(c+1) - 1),1);
    % Find node indices of the faces of the cell subset through indrection
    % map.
    ni = mcolon(G.faces.nodePos(cf), ...
                G.faces.nodePos(cf+1)-1)';
    
    W = [rldecode(cellno, nnode(cf)), G.faces.nodes(ni)];
    W = rlencode(sortrows(W));
    pos = cumsum([1; accumarray(W(:,1),1)]);
    n = W(:,2);

    assert (numel(pos) == numel(c) + 1);
end
