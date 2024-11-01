function varargout = invv(varargin)
% INVV(A, sz)
% computes inverse of sz(k)-by-sz(k) systems using lapack.
%
% The output is algorithmically equivalent to
%
%  p = cumsum([1, sz.^2]);
%  B = zeros([p(end) - 1, 1]);
%  for i = 1 : numel(sz),
%     ix    = p(i) : p(i + 1) - 1;
%     B(ix) = inv(reshape(A(ix), sz([i, i])));
%  end
%
% except for round-off errors.

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

if ispc,

   e = mexext('all');
   a = e(strcmp({ e.ext }, mexext)).arch;

   mwlib = @(lib) ...
      fullfile(matlabroot, 'extern', 'lib', a, 'microsoft', ...
               ['libmw', lib, '.lib']);

elseif isunix,

   mwlib = @(lib) ['-lmw', lib];

end

buildmex('invv.c', mwlib('lapack'), mwlib('blas'));

% Now, run it.
[varargout{1:nargout}] = invv(varargin{:});
