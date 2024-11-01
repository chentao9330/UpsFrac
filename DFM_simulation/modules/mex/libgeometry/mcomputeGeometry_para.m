function G = mcomputeGeometry(G,flag)
%Compute geometric primitves using compiled C code.
%
% SYNOPSIS:
%   G = mcomputeGeometry(G)
%
% PARAMETERS:
%   G - A grid_structure.
%
% RETURNS:
%   G - An updated grid_structure containing areas, volumes, normals and
%       centroids.
%
% SEE ALSO:
%   computeGeometry, grid_structure.

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
if(nargin==1)
  flag=0;
end
   if numel(G) > 0,
      if ~isfield(G(1), 'type'),
         warning(msgid('GridType:Unknown'),                         ...
                ['Input grid has no known type. ',                  ...
                 'I''ll assume it arose from the primordial soup...']);

         [ G(:).type ] = deal( {'Primordial Soup'} );
      end

      for k = 1 : numel(G),
         if size(G(k).nodes.coords, 2) ~= 3,
            error(['Function ''%s'' is only supported in three ', ...
                   'space dimensions.'], mfilename);
         end
         if(flag==0)
             [fa, fc, fn, cc, cv] = mex_compute_geometry(G(k));
          else
              [fa, fc, fn, cc, cv] = mex_compute_geometry_para(G(k),flag);   
         end

         G(k).faces.areas     = fa;
         G(k).faces.centroids = fc .';
         G(k).faces.normals   = fn .';
         G(k).cells.centroids = cc .';
         G(k).cells.volumes   = cv;

         G(k).type = [ G(k).type, { mfilename } ];
      end
   end
end
