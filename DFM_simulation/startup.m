function startup
%Amend MATLAB PATH to handle MRST implementation.

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

   d = fileparts(mfilename('fullpath'));
   m = fullfile(d, 'modules');

   p = split_path(genpath(d));

   i =     strncmp(m, p, length(m));
   i = i | ~cellfun(@isempty, regexp(p, '\.(git|hg|svn)'));
   i = i | cellfun(@isempty, p);

   addpath(p{~i});

   % Activate selected modules
   mrstModule add mex
end

%--------------------------------------------------------------------------

function p = split_path(p)
   try
      p = regexp(p, pathsep, 'split');
   catch  %#ok
      % Octave compatiblity.  It is an error to get here in an M run.
      p = strsplit(p, pathsep);
   end
end
