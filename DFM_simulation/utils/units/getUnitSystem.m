function u = getUnitSystem(s)
%Define unit conversion factors for input data.
%
% SYNOPSIS:
%   u = getUnitSystem(s)
%
% PARAMETERS:
%   s - String naming particular ECLIPSE deck unit conventions.  Accepted
%       values are:
%          - 'METRIC' --  Metric unit conventions (lengths in meter &c)
%          - 'FIELD'  --  Field unit conventions (lengths in feet &c)
%          - 'SI'     --  Strict SI conventions (all conversion factors 1)
%
% RETURNS:
%   u - Data structure of unit conversion factors from input data to MRST's
%       strict SI conventions.  Currently defines the following fields:
%
%           length    - Input unit of measurement for lengths/depths &c
%           time      - Input unit of measurement for time
%           viscosity - Input unit of measurement for viscosities
%           perm      - Input unit of measurement for permeabilities
%           trans     - Input unit of measurement for transmissibilitities
%
% SEE ALSO:
%   readGRDECL, convertInputUnits, convertFrom.

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

   switch lower(s),
      case 'metric'
         u = struct('length'   , meter             , ...
                    'time'     , day               , ...
                    'viscosity', centi*poise       , ...
                    'perm'     , milli*darcy       , ...
                    'trans'    , centi*poise * meter^3 / (day * barsa));

      case 'field',
         u = struct('length'   , ft                , ...
                    'time'     , day               , ...
                    'viscosity', centi*poise       , ...
                    'perm'     , milli*darcy       , ...
                    'trans'    , centi*poise * stb / (day * psia));

      case 'si',
         % SI units.  MRST extension.  Idempotency.
         u = struct('length'   , 1, ...
                    'time'     , 1, ...
                    'viscosity', 1, ...
                    'perm'     , 1, ...
                    'trans'    , 1);

      otherwise,
         error('Input unit system must be either METRIC, FIELD or SI.');
   end
end
