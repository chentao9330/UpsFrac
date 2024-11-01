function plotVelocityVE(g, state)
%Plot two-dimensional velocity vectors in solution.
%
% SYNOPSIS:
%   plotVelocityVE(G, state)
%
% PARAMETERS:
%   G     - A grid structure.
%
%   state - An MRST solution structure containing a valid flux field.
%
% RETURNS:
%   Nothing, but displays a QUIVER plot of the lateral velocity components.
%
% SEE ALSO:
%   quiver, grid_structure.

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

cellflux     = faceFlux2cellFlux(g, state.flux);

cellvelocity = cellFlux2cellVelocity(g, cellflux);

cellvelocity = bsxfun(@rdivide, cellvelocity, ...
                      sqrt(sum(cellvelocity .^ 2, 2)));

quiver(g.cells.centroids(:,1), ...
       g.cells.centroids(:,2), ...
       cellvelocity(:,1), cellvelocity(:,2));
end
%{
return
X=reshape(g.cells.centroids(:,1),g.cartDims);
Y=reshape(g.cells.centroids(:,2),g.cartDims);
vx=reshape(cellvel(:,1),g_top.cartDims);vy=reshape(cellvel(:,2),g_top.cartDims);
v_norm=sqrt(vx.^2+vy.^2);
pos=[g_top.cells.centroids(startcells,:)]%positions to start streamlines
r=100;vv=0:01:pi;clf,streamline(X',Y',(vx./v_norm)',(vy./v_norm)',pos(:,1),pos(:,2))
%}
