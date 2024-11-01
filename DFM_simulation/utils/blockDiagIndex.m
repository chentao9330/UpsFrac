function varargout = blockDiagIndex(m, n)
%Compute subscript or linear index to nonzeros of block-diagonal matrix
%
% SYNOPSIS:
%   [i, j] = blockDiagIndex(m,n)
%   [i, j] = blockDiagIndex(m)
%
% REQUIRED PARAMETERS:
%   m, n,  - vectors of dimensions of each diagonal block. 
%            
% RETURNS:
%   i,j    - index vectors that may be used to form sparse block-diagonal
%            matrix.  See example.
%
% EXAMPLE
%
%   m     = [1;2;3]; 
%   n     = [2;3;4];    
%   [i,j] = blockDiagIndex(m,n);
%   A     = sparse(i,j,1:sum(m.*n));
%   full(A)
%
% SEE ALSO:
%   rldecode, mcolon.

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

   if nargin == 1, n = m; end
   pos = cumsum([1;double(m(:))]);
   p1  = pos(1:end-1);
   p2  = pos(2:end)-1;
   i   = mcolon(rldecode(p1, n(:)),  rldecode(p2, n(:)))';
  
   pos = cumsum([1;double(n(:))]);   
   p1  = pos(1:end-1);
   p2  = pos(2:end)-1;
   j   = rldecode(mcolon(p1, p2)', rldecode(m(:),n(:)));
   
   
   if nargout < 2,
      varargout{1} = sub2ind([sum(m), sum(n)], i, j);
   elseif nargout == 2,
      varargout{1} = i;
      varargout{2} = j;
   else
      error('Huh!? Too many output arguments');
   end
end
