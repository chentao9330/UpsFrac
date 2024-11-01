function [kro, Soco, Socr, Somax] = sof2(T, varargin)
% Construct oil-water or oil-gas relperm eval. functions from SOF2 table.
%
% SYNOPSIS:
%   [kro, Soco, Socr, Somax] = sof3(table)
%
% PARAMETERS:
%   table - Saturation function table as defined through the input deck
%           keyword 'SOF2'.  May be a numeric array or a a cell array of
%           such arrays.
%
% RETURNS:
%   In the following, 'ntab' refers to the number of input tables
%   represented by the input parameter 'table'.  If 'table' is a single
%   numeric array, then ntab = 1, otherwise (i.e., when 'table' is a cell
%   array of numeric arrays) ntab = numel(table).
%
%   In an ECLIPSE input deck, 'ntab' typically corresponds to the 'NTSFUN'
%   parameter specified in item one of the 'TABDIMS' keyword.
%
%   kro    - 
% 
%   Soco   - 
% 
%   Socr   - 
% 
%   Somax  - 
% 
%
% SEE ALSO:
%   readRelPermTable, swfn, sgfn, sof2, swof, sgof, initBlackOilRelPerm.

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
   assert (isnumeric(T) || (iscell(T) && all(cellfun(@isnumeric, T))),  ...
          ['Parameter ''table'' must be a numeric array or cell array', ...
           ' of numeric arrays.']);

   if isnumeric(T), T = { T }; end

   ntab = numel(T);

   kro  = cell ([1, ntab]);
   Soco  = zeros([1, ntab]);
   Socr  = zeros([1, ntab]);
   Somax = zeros([1, ntab]);
   for k = 1 : ntab,
      [Soco(k), Socr(k), Somax(k)] = table_points(T{k});

      if Soco(k) > 0, T{k} = [0, T{k}(1, 2:end); T{k}]; end

      kro{k} = @(so) interp_oil(so, T{k}, Soco(k));
   end
end
      
%--------------------------------------------------------------------------
function [Soco, Socr, Somax] = table_points(T)
   % Basic validation of input data.
   assert (~any(any(T(:, 1:2) < 0)), ...
      'Saturation and rel-perm values must be non-negative.');
   

   assert (~ (T( 1 ,2) > 0));        % krow(Soco)     == 0

   assert (all (diff(T(:,1)) > 0));  % Sw monotonically increasing
   assert (~any(diff(T(:,2)) < 0));  % Level or increasing down column   
   
   %-----------------------------------------------------------------------
   % Data OK.  Extract critical values from table.
   %
   
   % 1) Connate water saturation (>= 0) is first Sw encountered in table.
   %
   Soco = T(1,1);  assert (~ (Soco < 0));
   
   % 2) Residual oil saturation in oil-water system
   i    = find(T(:,2) > 0, 1, 'first');  assert (~isempty(i) && (i > 1));
   Socr = T(i - 1, 1);                   assert (~ (Socr < Soco));

   % 3) Maximal water saturation in table
   Somax = max(T(:,1));
end

%--------------------------------------------------------------------------
% Water relative permeability as a function of water saturation, sw.
% Note: Water saturation is always >= Swco.
function varargout = interp_oil(so, T, Soco)
   so = max(so, Soco);
   
   varargout{1}    = interpTable(T(:,1), T(:,2), so);
   
   if nargout > 1,
      varargout{2} = dinterpTable(T(:,1), T(:,2), so);
   end
end
