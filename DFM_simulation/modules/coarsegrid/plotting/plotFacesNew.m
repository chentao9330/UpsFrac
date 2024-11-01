function h = plotFacesNew(G, varargin)
  %Plot selection of coloured grid faces to current axes (reversed Z axis).
%
% SYNOPSIS:
%       plotFaces(G, faces)
%       plotFaces(G, faces, 'pn1', pv1, ...)
%       plotFaces(G, faces, colour)
%       plotFaces(G, faces, colour, 'pn1', pv1, ...)
%   h = plotFaces(...)
%
% PARAMETERS:
%   G       - Grid data structure.
%
%   faces   - Vector of face indices.  The graphical output of 'plotFaces'
%             will be restricted to the subset of grid faces from 'G'
%             represented by 'faces'.
%
%   colour  - Colour data specification.  Either a MATLAB 'ColorSpec'
%             (i.e., an RGB triplet (1-by-3 row vector) or a short or long
%             colour name such as 'r' or 'cyan'), or a PATCH
%             'FaceVertexCData' table suiteable for either indexed or
%             'true-colour' face colouring.  This data *MUST* be an m-by-1
%             column vector or an m-by-3 matrix.  We assume the following
%             conventions for the size of the colour data:
%
%                - ANY(SIZE(colour,1) == [1, NUMEL(faces)])
%                  One (constant) indexed colour for each face in 'faces'.
%                  This option supports 'flat' face shading only.  If
%                  SIZE(colour,1) == 1, then the same colour is used for
%                  all faces in 'faces'.
%
%                - SIZE(colour,1) == G.nodes.num
%                  One (constant) indexed colour for each node in 'faces'.
%                  This option must be chosen in order to support
%                  interpolated face shading.
%
%             OPTIONAL.  Default value: colour = 'y' (shading flat).
%
%   'pn'/pv - List of other property name/value pairs.  OPTIONAL.
%             This list will be passed directly on to function PATCH
%             meaning all properties supported by PATCH are valid.
%
% RETURNS:
%   h - Handle to resulting PATCH object.  The patch object is added to the
%       current AXES object.
%
% NOTES:
%   Function 'plotFaces' is implemented directly in terms of the low-level
%   function PATCH.  If a separate axes is needed for the graphical output,
%   callers should employ function newplot prior to calling 'plotFaces'.
%
% EXAMPLE:
%   % Plot grid with boundary faces on left side in red colour:
%   G     = cartGrid([5, 5, 2]);
%   faces = boundaryFaceIndices(G, 'LEFT', 1:5, 1:2, []);
%   plotGrid (G, 'faceColor', 'none'); view(3)
%   plotFaces(G, faces, 'r');
%
% SEE ALSO:
%   plotCellData, plotGrid, newplot, patch, shading.


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

   %% Split input and options, assuming first string input is an option
   i           = cellfun(@ischar, varargin);
   first       = find(i, 1, 'first');
   if mod(numel(varargin)-first+1, 2), first = first + 1; end
   if any(first),
      other_input = varargin(1:first-1);
      varargin    = varargin(first:end);
   else
      other_input = varargin;
      varargin = {};
   end   
%   if mod(numel(varargin), 2) && numel(varargin{1}) == 1, 
%      varargin = {'edgecolor', varargin{1}, 'facecolor', varargin{:}};
%   end
   
   % There should be an even number of elements in varargin here:
   assert(~mod(numel(varargin), 2), 'Huh?!');
   
   color = 'y';
   switch numel(other_input),
      case 0
         if G.griddim == 2, 
            faces = 1:G.faces.num;
         else
            faces = boundaryFaces(G, 1:G.cells.num);
         end
         if ~any(strcmpi(varargin,'FaceColor')),
            varargin = [varargin {'FaceColor','y'}];
         end
      case 1
         faces = other_input{1};
         if ~any(strcmpi(varargin,'FaceColor')), 
            varargin = [varargin {'FaceColor','y'}];
         end
      case 2
         faces = other_input{1};
         color = other_input{2};
      otherwise
         error('What!?');
   end
   
   
   assert(min(faces) > 0,'Cannot plot zero or negativ face numbers');
   assert(max(faces) <= G.faces.num,...
      'Faces in face list exceed number of faces in grid.');

   if G.griddim == 3,
      % If G is a coarsegrid, lookup finegrid data in parent
      if isfield(G, 'parent'),
         [f, fno] = getSubFaces(G, faces);
         if mod(numel(varargin), 2), varargin{1} = varargin{1}(fno);end
        
         %% Marker-related otions are collected in marker_opts
         ix          = rldecode(strncmpi('marker', varargin(1:2:end), 6), 2);
         marker_opts = varargin(ix);
         varargin    = varargin(~ix);


         %% Edge-related options are collected in edge_opts 
         ix = rldecode(strncmpi('edge', varargin(1:2:end), 4)', 2) | ...
            rldecode(strncmpi('line', varargin(1:2:end), 4)', 2);
         edge_opts = varargin(ix);
         varargin  = varargin(~ix);
         
         
         h = plotPatches(G.parent, f, 'edgec', 'none', varargin{:});
         set(get(h, 'Parent'), 'ZDir', 'reverse')
         h = [h; plotFaceOutline(G, faces, edge_opts{:})];
         
         if numel(marker_opts) > 0,
            if isfield(G, 'parent'),
               cg = G;
               % [f, fno]  = getSubFaces(G, faces);
               G  = G.parent;
            else
               f  = faces;
            end
            
            % try
            %   faceno = rldecode(faces(:), (cg.faces.connPos(faces+1)-cg.faces.connPos(faces))*2, 1);               
            % catch
            %   faceno = rldecode(faces(:), 2);
            % end
            
            ff = cg.faces.fconn;
            ffno = rldecode(1:cg.faces.num, diff(cg.faces.connPos), 2)';
            
            % Bug here:  A node is shared by two (or more) faces in 2D and
            % three (or more faces) in 3D.  This definition does not
            % coincide with "shared by two edges" hack above.
            
            % This check must be done per block: Check implementation of
            % cellNodes!
            
            [d, p] = copyRowsFromPackedData(G.faces.nodes, G.faces.nodePos, ff);
            fedges = [d, d(rot(p, 1))];
            faceno = rldecode(ffno, diff(p));
            tmp    = rlencode(sortrows([faceno, fedges(:,1); faceno, fedges(:,2)], [2,1]));
            [n,n]  = rlencode(tmp(:,2)); N=rldecode(n,n);
            nodes  = rlencode(tmp(N>2, :));
            
            holdstate   = ishold;
            hold on;
            %plot3(G.nodes.coords(nodes,1), G.nodes.coords(nodes, 2), G.nodes.coords(nodes, 3), 'linestyle','none', marker_opts{:});
            if ~holdstate, hold off; end         
            
            %warning('Nodes in 3d is currently not supported');
         end

      else
         h = plotPatches(G, faces, color, 'EdgeColor', 'k', varargin{:});
         set(get(h, 'Parent'), 'ZDir', 'reverse')
      end 
   else
      % If G is a coarsegrid, lookup finegrid data in parent
      if isfield(G, 'parent'),
         cg = G;
         f  = getSubFaces(G, faces);
         G  = G.parent;
      else
         f  = faces;
      end


      %% Separate otions: marker-related stuff is sent to separate plotting
      ix          = rldecode(strncmpi('marker', varargin(1:2:end), 6), 2);
      marker_opts = varargin(ix);
      varargin    = varargin(~ix);
      
      %% Plot fine grid edges specified by either coarse or fine grid.
      ix          = mcolon(G.faces.nodePos(f), G.faces.nodePos(f+1)-1);
      edges       = reshape(G.faces.nodes(ix), 2, [])';
      h           = plotLineSegments(G, edges, varargin{:});

      %% Find unique endpoints - if varargin contains 'marker*'      
      if numel(marker_opts) > 0
         try
            faceno = rldecode(faces(:), (cg.faces.connPos(faces+1)-cg.faces.connPos(faces))*2, 1);
         catch %#ok
            faceno = rldecode(faces(:), 2);
         end
         e     = reshape(G.faces.nodes, 2, [])';
         e     = reshape(e(f,:)', [], 1);
         [E,n] = rlencode(sortrows([faceno, e]));
         nodes = E(n==1,2);
         
         holdstate   = ishold;
         hold on;
         patch('vertices',G.nodes.coords(nodes,:), 'faces',(1:numel(nodes))',marker_opts{:});
         if ~holdstate, hold off; end
      end
   end

end

function h = plotLineSegments(G, e, varargin)
% Plot all line segments given by node pairs in each row in e.
   e = unique([e;fliplr(e)], 'rows');
   h = patch('vertices', G.nodes.coords, 'faces', e, varargin{:});
end
function [subf, fno] = getSubFaces(G, f)
   ix   = mcolon(G.faces.connPos(f), G.faces.connPos(f+1)-1);
   subf = G.faces.fconn(ix);
   fno  = rldecode(1:numel(f), G.faces.connPos(f+1)-G.faces.connPos(f), 2)';
end
function ix = rot(pos, offset)
   num    = diff(pos);
   offset = mod(offset, num); % net offset
   ix     = zeros(max(pos)-1, 1);

   ix(mcolon(pos(1:end-1), pos(1:end-1)+num-offset-1)) = ...
      mcolon(pos(1:end-1)+offset, pos(2:end)-1);
   ix(mcolon(pos(1:end-1)+num-offset, pos(2:end)-1)) = ...
      mcolon(pos(1:end-1), pos(2:end)-1-num+offset);
end
function [d, p] = copyRowsFromPackedData(d, p, rows)
   d = d(mcolon(p(rows), p(rows+1)-1));
   p = cumsum([1;double(p(rows+1)-p(rows))]);
end
   
