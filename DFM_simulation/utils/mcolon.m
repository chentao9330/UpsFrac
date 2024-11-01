function x = mcolon(lo, hi, s)
%Compute vector of consecutive indices from separate lower/upper bounds.
%
% SYNOPSIS:
%   ind = mcolon(lo, hi)
%   ind = mcolon(lo, hi, stride)
%
% PARAMETERS:
%   lo  - Vector of start values (lower bounds).
%   hi  - Vector of end values (upper bounds).
%   s   - Vector of strides.
%         Optional.  Default value: s = ones([numel(lo), 1]) (unit stride).
%
% RETURNS:
%   ind - [lo(1):hi(1)     , lo(2):hi(2)     , ..., lo(end):hi(end)       ]
%   ind - [lo(1):s(1):hi(1), lo(2):s(2):hi(2), ..., lo(end):s(end):hi(end)]
%
%   Note that 'ind' has type DOUBLE irrespective of the type of its input
%   parameters.
%
% EXAMPLE:
%   lo  = [1 1 1 1]; hi = [2 3 4 5];
%   ind = mcolon(lo, hi)
%
% NOTE:
%   MCOLON may be implemented in terms of ARRAYFUN and HORZCAT, e.g.,
%      ind = arrayfun(@colon, lo, hi, 'UniformOutput', false);
%      ind = [ind{:}];
%   or
%      ind = arrayfun(@colon, lo, s, hi, 'UniformOutput', false);
%      ind = [ind{:}];
%   but the current implementation is faster.

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

% Support general input class.
lo = double(lo);
hi = double(hi);

if numel(lo) == 1,
   lo = repmat(lo, size(hi));
end
if numel(hi)==1,
   hi = repmat(hi, size(lo));
end   

if numel(lo) ~= numel(hi)
  error('In mcolon: lo and hi must have same number of elements!');
elseif numel(lo) == 0,
  x = [];
elseif nargin < 3
   
  % Remove lo-hi pairs where numel(lo:hi)==0
  i    = hi>=lo;
  hi   = hi(i);
  lo   = lo(i);
  if sum(i) == 0, x=[];return; end
  m    = numel(lo);
  d    = double(hi - lo + 1);
  n    = sum(d);

  x    = ones(1,n);
  x(1) = lo(1);
  x(1+cumsum(d(1:end-1))) = lo(2:m) - hi(1:m-1);
  x    = cumsum(x);
else
  s = double(s);
  if numel(s) == 1,
     s = repmat(s, size(lo));
  end
  
  % Remove lo-hi-s triplets where numel(lo:s:hi)==0
  i    = ((hi >= lo) & (s > 0)) | ((hi <= lo) & (s < 0));
  hi   = hi(i);
  lo   = lo(i);
  s    = s(i);

  if sum(i) == 0, x=[];return; end

  % Compute lo + (0:(hi-lo)/stride)*stride
  % Fix or hack: avoid roundoff error in floor((hi-lo)./s) when hi-lo = N*s
  % for some natural number N.
  e    =  (1- 2*(hi<lo)) * eps;
  hi   = fix((e+hi-lo)./s);


  m    = numel(lo);
  d    = double(hi + 1);
  n    = sum(d);

  assert (all(d > 0), ...
         ['Internal error in ''%s'': Bins with non-positive ', ...
          'number of elements detected'], mfilename);

  ind = 1+cumsum(d(1:end-1));

  % Expand lo  to [lo(1) lo(1) ... lo(2) lo(2) ... lo(end)]
  LO   = zeros(1,n);
  LO(1)=lo(1);
  LO(ind) = lo(2:m) - lo(1:m-1);
  LO   = cumsum(LO);

  % Expand stride
  S    = zeros(1,n);
  S(1) = s(1);
  S(ind) = s(2:m) - s(1:m-1);
  S    = cumsum(S);

  x    = ones(1,n);
  x(1) = 0;
  x(ind) = -hi(1:m-1);

  x    = cumsum(x).*S + LO;
end
end
