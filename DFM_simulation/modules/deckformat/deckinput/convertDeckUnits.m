function deck = convertDeckUnits(deck, varargin)
%Convert ECLIPSE/FrontSim input deck units to MRST conventions
%
% SYNOPSIS:
%   deck = convertDeckUnits(deck)
%   deck = convertDeckUnits(deck, 'pn1', pv1, ...)
%
% PARAMETERS:
%   deck    - An ECLIPSE/FrontSim input deck as defined by function
%             'readEclipseDeck'.
%
%   'pn'/pv - A list of 'key'/value pairs defining optional parameters.
%             The supported options are:
%               - verbose -- Whether or not to emit informational messages
%                            if a section is skipped.
%                            Logical.  Default value: verbose = mrstVerbose
%
% RETURNS:
%   deck - An ECLIPSE/FrontSim deck where all quantities affected by
%          varying units of measurements have been converted to MRST's
%          strict, SI only, conventions.
%
% SEE ALSO:
%   readEclipseDeck, mrstVerbose.

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

   opt = struct('verbose', mrstVerbose);
   opt = merge_options(opt, varargin{:});

   assert (isstruct(deck) && isfield(deck, 'RUNSPEC'), ...
           'Input does not appear to be a valid deck');

   u = unit_system(deck.RUNSPEC);

   for sect = reshape(fieldnames(deck), 1, []),
      s = sect{1};

      switch s,
         case {'RUNSPEC', 'GRID', 'PROPS', 'SOLUTION', 'SCHEDULE'},
            cvrt     = str2func(['convert', s]);
            deck.(s) = cvrt(deck.(s), u);

         case 'UnhandledKeywords',
            continue   % MRST specific

         otherwise
            dispif(opt.verbose, ...
                   'No converter needed in section ''%s''.\n', s);
      end
   end
end

%--------------------------------------------------------------------------
% Private helpers follow
%--------------------------------------------------------------------------

function u = unit_system(rspec)
   metric = isfield(rspec, 'METRIC');
   field  = isfield(rspec, 'FIELD');
   lab    = isfield(rspec, 'LAB');
   SI     = isfield(rspec, 'SI');

   if sum([metric, field, lab, SI]) ~= 1,
      error(id('USys:Unknown'), ...
            'Input unit system must be either METRIC, FIELD, LAB, or SI.');
   end

   if metric,
      u = struct('length'   , meter             , ...
                 'time'     , day               , ...
                 'density'  , kilogram / meter^3, ...
                 'press'    , barsa             , ...
                 'concentr' , kilogram / meter^3, ... % Concentration
                 'compr'    , 1 / barsa         , ... % Compressibility
                 'viscosity', centi*poise       , ...
                 'perm'     , milli*darcy       , ...
                 'liqvol_s' , meter^3           , ... % Liquid vol, surf
                 'liqvol_r' , meter^3           , ... % Liquid vol, res
                 'gasvol_s' , meter^3           , ... % Gas vol, surf
                 'gasvol_r' , meter^3           , ... % Gas vol, res
                 'trans'    , centi*poise * meter^3 / (day * barsa));
   elseif field
      u = struct('length'   , ft                , ...
                 'time'     , day               , ...
                 'density'  , pound / ft^3      , ...
                 'press'    , psia              , ...
                 'concentr' , pound / stb       , ... % Concentration
                 'compr'    , 1 / psia          , ...
                 'viscosity', centi*poise       , ...
                 'perm'     , milli*darcy       , ...
                 'liqvol_s' , stb               , ...
                 'liqvol_r' , stb               , ...
                 'gasvol_s' , 1000 * ft^3       , ... % Mscf
                 'gasvol_r' , stb               , ...
                 'trans'    , centi*poise * stb / (day * psia));
   elseif lab,
      u = struct('length'   , centi*meter           , ...
                 'time'     , hour                  , ...
                 'density'  , gram / (centi*meter)^3, ...
                 'press'    , atm                   , ...
                 'concentr' , gram / (centi*meter)^3, ...
                 'compr'    , 1 / atm               , ...
                 'viscosity', centi*poise           , ...
                 'perm'     , milli*darcy           , ...
                 'liqvol_s' , (centi*meter)^3       , ...
                 'liqvol_r' , (centi*meter)^3       , ...
                 'gasvol_s' , (centi*meter)^3       , ...
                 'gasvol_r' , (centi*meter)^3       , ...
                 'trans'    , centi*poise * (centi*meter)^3 / (hour * atm));
   else
      % SI units.  MRST extension.  Idempotency.
      u = struct('length'   , 1, ...
                 'time'     , 1, ...
                 'density'  , 1, ...
                 'press'    , 1, ...
                 'concentr' , 1, ...
                 'compr'    , 1, ...
                 'viscosity', 1, ...
                 'perm'     , 1, ...
                 'liqvol_s' , 1, ...
                 'liqvol_r' , 1, ...
                 'gasvol_s' , 1, ...
                 'gasvol_r' , 1, ...
                 'trans'    , 1);
   end
end

%--------------------------------------------------------------------------

function rspec = convertRUNSPEC(rspec, varargin)                       %#ok
   usys = { 'METRIC', 'FIELD', 'LAB' };

   rspec = rmfield(rspec, usys(isfield(rspec, usys)));

   rspec.SI = true;
end

%--------------------------------------------------------------------------

function grid = convertGRID(grid, u)                                   %#ok
   for kw = reshape(fieldnames(grid), 1, []),
      key = kw{1};

      switch key,
         case {'PERMX' , 'PERMXY', 'PERMXZ', ...
               'PERMYX', 'PERMY' , 'PERMYZ', ...
               'PERMZX', 'PERMZY', 'PERMZ' , ...
               'PERMXX', 'PERMYY', 'PERMZZ'},
            grid.(key) = convertFrom(grid.(key), u.perm);

         case {'DXV'   , 'DYV'   , 'DZV'   , 'DEPTHZ', ...
               'DX'    , 'DY'    , 'DZ'    , 'TOPS'  , ...
               'COORDX', 'COORDY', 'COORDZ',           ...
               'COORD' , 'ZCORN'                     },
            grid.(key) = convertFrom(grid.(key), u.length);

         case 'MAPAXES',
            unt = u.length;
            if isfield(grid, 'MAPUNITS'),
               switch grid.MAPUNITS,
                  case 'METRES', unt = meter;
                  case 'FEET'  , unt = ft;
                  otherwise
                     warning('MapUnits:Unknown', ...
                             'Unknown map units ''%s''. Using METRES', ...
                             grid.MAPUNITS);
               end
            end
            grid.(key) = convertFrom(grid.(key), unt);

         case 'MAPUNITS',
            continue;  % Handled in 'MAPAXES'.

         case {'MINPV', 'MINPVV'},
            grid.(key) = convertFrom(grid.(key), u.liqvol_r);

         case 'NNC',
            grid.(key)(:,7) = convertFrom(grid.(key)(:,7), u.trans);

         case {'TRANX', 'TRANY', 'TRANZ'},
            grid.(key) = convertFrom(grid.(key), u.trans);

         case 'PINCH',
            i    = [1, 3];
            data = convertFrom([ grid.(key){i} ], u.length);

            grid.(key)(i) = num2cell(data);   clear i data

         case {'cartDims',                        ...  % MRST specific
               'ACTNUM', 'NTG', 'PORO', 'MULTPV', ...
               'FLUXNUM',                         ...
               'FAULTS', 'MULTFLT',               ...
               'MULTX' , 'MULTX_', 'MULTY', 'MULTY_', 'MULTZ', 'MULTZ_'},
            continue;  % Dimensionless

         otherwise
            error(id('GRID:NoConverter'), ...
                  'No known unit conversion in GRID for ''%s''.', key);
      end
   end
end

%--------------------------------------------------------------------------

function props = convertPROPS(props, u)                                %#ok
   for kw = reshape(fieldnames(props), 1, []),
      key = kw{1};

      switch key,
         case 'DENSITY',
            props.(key) = convertFrom(props.(key), u.density);

         case 'RKTRMDIR'
            continue; %dimensionless

         case 'ROCK',
            unt         = [u.press, u.compr, u.compr, u.compr, 1, 1];
            unt         = unt(1 : size(props.(key), 2));
            props.(key) = convertFrom(props.(key), unt);

         case 'ROCKTAB',
            unt         = [u.press, ones(1, size(props.(key){1}, 2)-1)];
            for t = 1:numel(props.(key))
                props.(key){t} = convertFrom(props.(key){t}, unt);
            end

         case 'PLYMAX',
            unt         = [u.concentr, u.concentr];
            props.(key) = convertFrom(props.(key), unt);

         case 'PLYROCK',
            unt         = [1, 1, u.concentr, 1, 1];
            props.(key) = convertFrom(props.(key), unt);

         case {'PLYVISC', 'PLYADS'},
            unt         = [u.concentr, 1];
            for t = 1 : numel(props.(key)),
                props.(key){t} = convertFrom(props.(key){t}, unt);
            end

         case 'PVCDO',
            unt         = [u.press, 1, u.compr, u.viscosity, u.compr];
            props.(key) = convertFrom(props.(key), unt);

         case 'PVDG',
            unt = [u.press, u.gasvol_r/u.gasvol_s, u.viscosity];
            for t = 1 : numel(props.(key)),
               props.(key){t} = convertFrom(props.(key){t}, unt);
            end

         case 'PVDO',
            assert (iscell(props.(key)));

            unt = [u.press, u.liqvol_r/u.liqvol_s, u.viscosity];
            for t = 1 : numel(props.(key)),
               props.(key){t} = convertFrom(props.(key){t}, unt);
            end

         case 'PVTG',
            assert (iscell(props.(key)));

            uk = u.press;
            ud = [u.liqvol_s/u.gasvol_s, ...
                  u.gasvol_r/u.gasvol_s, u.viscosity];

            for t = 1 : numel(props.(key)),
               props.(key){t}.key  = convertFrom(props.(key){t}.key , uk);
               props.(key){t}.data = convertFrom(props.(key){t}.data, ud);
            end

         case 'PVTO',
            uk = u.gasvol_s / u.liqvol_s;
            ud = [u.press, u.liqvol_r/u.liqvol_s, u.viscosity];

            for t = 1 : numel(props.(key)),
               props.(key){t}.key  = convertFrom(props.(key){t}.key , uk);
               props.(key){t}.data = convertFrom(props.(key){t}.data, ud);
            end

         case 'PVTW',
            unt         = [u.press, u.liqvol_r/u.liqvol_s, ...
                           u.compr, u.viscosity, u.compr];
            props.(key) = convertFrom(props.(key), unt);

         case 'RSCONSTT'
            unt         = [u.gasvol_s/u.liqvol_s, u.press];
            props.(key) = convertFrom(props.(key), unt);

         case {'SGFN', 'SWFN'},
            unt = [1, 1, u.press];

            for t = 1 : numel(props.(key)),
               props.(key){t} = convertFrom(props.(key){t}, unt);
            end

         case {'SGOF', 'SWOF', 'SGWFN', 'SLGOF'},
            assert (iscell(props.(key)));

            unt = [1, 1, 1, u.press];

            for t = 1 : numel(props.(key)),
               props.(key){t} = convertFrom(props.(key){t}, unt);
            end

         case {'SOF2', 'SOF3', 'STONE', 'STONE1', 'STONE2', ...
               'SIMPLE', 'TLMIXPAR', ...
               ...
               'SWL'   ,            'ISWL' ,           ...
               'SWLX'  , 'SWLX-'  , 'ISWLX', 'ISWLX-', ...
               'SWLY'  , 'SWLY-'  , 'ISWLY', 'ISWLY-', ...
               'SWLZ'  , 'SWLZ-'  , 'ISWLZ', 'ISWLZ-', ...
               'SWCR'  ,            'SWU'  ,           ...
               'SWCRX' , 'SWCRX-' , 'SWUX' , 'SWUX-' , ...
               'SWCRY' , 'SWCRY-' , 'SWUY' , 'SWUY-' , ...
               'SWCRZ' , 'SWCRZ-' , 'SWUZ' , 'SWUZ-' , ...
               'ISWCR' ,            'ISWU' ,           ...
               'ISWCRX', 'ISWCRX-', 'ISWUX', 'ISWUX-', ...
               'ISWCRY', 'ISWCRY-', 'ISWUY', 'ISWUY-', ...
               'ISWCRZ', 'ISWCRZ-', 'ISWUZ', 'ISWUZ-'},
            continue;  % Dimensionless

         otherwise
            error(id('PROPS:NoConverter'), ...
                  'No known unit conversion in PROPS for ''%s''.', key);
      end
   end
end

%--------------------------------------------------------------------------

function soln = convertSOLUTION(soln, u)                               %#ok
   for kw = reshape(fieldnames(soln), 1, []),
      key = kw{1};

      switch key,
         case 'AQUCT',                               %  7  8  9 10 11
            unt = [1, u.length, u.press u.perm, 1,...
               1./u.press, u.length, u.length, 1,...
               1,1,0.0,NaN];
            soln.(key) = convertFrom(soln.(key), unt);

         case 'AQUANCON'
            soln.(key)(:,9) = ...
               cellfun(@(c) convertFrom(c, u.length ^ 2), ...
                       soln.(key)(:,9), 'UniformOutput', false);

         case 'EQUIL',                               %  7  8  9 10 11
            unt = [repmat([u.length, u.press], [1, 3]), 1, 1, 1, 1, 1];

            unt        = unt(1 : size(soln.(key), 2));
            soln.(key) = convertFrom(soln.(key), unt);

         case 'DATUM',
            soln.(key) = convertFrom(soln.(key), u.length);

         case {'PBVD', 'PDVD'},
            unt = [u.length, u.press];

            for reg = 1 : numel(soln.(key)),
               soln.(key){reg} = convertFrom(soln.(key){reg}, unt);
            end

         case {'PBUB', 'PRESSURE'}
               soln.(key) = convertFrom(soln.(key), u.press);

         case 'RSVD',
            unt = [u.length, u.gasvol_s/u.liqvol_s];

            for reg = 1 : numel(soln.(key)),
               soln.(key){reg} = convertFrom(soln.(key){reg}, unt);
            end

         case 'RVVD',
            unt = [u.length, u.liqvol_s/u.gasvol_s];

            for reg = 1 : numel(soln.(key)),
               soln.(key){reg} = convertFrom(soln.(key){reg}, unt);
            end

         case 'RS',
            soln.(key) = convertFrom(soln.(key), u.gasvol_s / u.liqvol_s);

         case 'RV',
            soln.(key) = convertFrom(soln.(key), u.liqvol_s / u.gasvol_s);

         case {'SGAS', 'SOIL', 'SWAT'},
            continue;  % Dimensionless

         otherwise
            error(id('SOLUTION:NoConverter'), ...
                  'No known unit conversion in SOLUTION for ''%s''.', key);
      end
   end
end

%--------------------------------------------------------------------------

function schd = convertSCHEDULE(schd, u)                               %#ok
   schd.step.val = convertFrom(schd.step.val, u.time);
   if(isfield(schd,'control')) 
    for c = 1 : numel(schd.control),
          schd.control(c) = convertControl(schd.control(c), u);
    end
   end
end

%--------------------------------------------------------------------------

function ctrl = convertControl(ctrl, u)
   for kw = reshape(fieldnames(ctrl), 1, []),
      key = kw{1};

      switch key,
         case 'COMPDAT',
            if ~isempty(ctrl.(key)),
               ctrl.(key) = convertCompDat(ctrl.(key), u);
            end

         case 'WCONHIST',
            if ~isempty(ctrl.(key)),
               ctrl.(key) = convertWconHist(ctrl.(key), u);
            end

         case 'WCONINJ',
            if ~isempty(ctrl.(key)),
               ctrl.(key) = convertWConInj(ctrl.(key), u);
            end

         case 'WCONINJE',
            if ~isempty(ctrl.(key)),
               ctrl.(key) = convertWConInje(ctrl.(key), u);
            end
            
         case 'WCONINJH',
            if ~isempty(ctrl.(key)),
               ctrl.(key) = convertWConInjh(ctrl.(key), u);
            end

         case 'WCONPROD',
            if ~isempty(ctrl.(key)),
               ctrl.(key) = convertWConProd(ctrl.(key), u);
            end
            
         case 'WPOLYMER',
            if ~isempty(ctrl.(key)),
               ctrl.(key) = convertWPolymer(ctrl.(key), u);
            end

         case 'GCONPROD',
            if ~isempty(ctrl.(key)),
               ctrl.(key) = convertGConProd(ctrl.(key), u);
            end

         case 'GCONINJE',
            if ~isempty(ctrl.(key)),
               ctrl.(key) = convertGconInje(ctrl.(key), u);
            end

         case 'WELSPECS',
            if ~isempty(ctrl.(key)),
               ctrl.(key) = convertWellSpecs(ctrl.(key), u);
            end

         case {'GRUPTREE', 'WGRUPCON'},
            continue; % Dimensionless

         case 'DRSDT',
            unt           = (u.gasvol_s / u.liqvol_s) / u.time;
            ctrl.(key){1} = convertFrom(ctrl.(key){1}, unt);

         otherwise
            error(id('SCHEDULE:NoConverter'), ...
                  'No known unit conversion in SCHEDULE for ''%s''.', key);
      end
   end
end

%--------------------------------------------------------------------------

function compdat = convertCompDat(compdat, u)
   c = [8      , 9       , 10             , 14      ];
   u = [u.trans, u.length, u.perm*u.length, u.length];

   for n = 1 : numel(c),
      compdat(:, c(n)) = cellfun(@(x) convertFrom(x, u(n)), ...
                                 compdat(:, c(n)), 'UniformOutput', false);
   end
end

%--------------------------------------------------------------------------

function wch = convertWconHist(wch, u)
   lrat = u.liqvol_s / u.time;
   grat = u.gasvol_s / u.time;

   c = [4   , 5   , 6   , 9      , 10     , 11  ];
   u = [lrat, lrat, grat, u.press, u.press, grat];

   c = c(c <= size(wch, 2));

   for n = 1 : numel(c),
      wch(:, c(n)) = cellfun(@(x) convertFrom(x, u(n)), ...
                             wch(:, c(n)), 'UniformOutput', false);
   end
end

%--------------------------------------------------------------------------

function wconinj = convertWConInj(wconinj, u)
   % Surface flow rate
   is_gas = strcmpi(wconinj(:, 2), 'gas');
   is_oil = strcmpi(wconinj(:, 2), 'oil');

   r_unit         = repmat(u.liqvol_s / u.time, [numel(is_gas), 1]);
   r_unit(is_gas) = u.gasvol_s / u.time;

   wconinj(:, 5)  = arrayfun(@convertFrom,                   ...
                             vertcat(wconinj{:, 5}), r_unit, ...
                             'UniformOutput', false);

   % Vapourised gas concentration
   if any(is_gas | is_oil),
      R_unit = repmat(u.liqvol_s / u.gasvol_s, [numel(is_gas), 1]);
      R_unit(is_oil) = u.gasvol_s / u.liqvol_s;

      wconinj(:, 12) = arrayfun(@convertFrom,                    ...
                                vertcat(wconinj{:, 12}), R_unit, ...
                                'UniformOutput', false);
   end

   % Handle default pressure
   i = vertcat(wconinj{:,9}) == -1;
   wconinj(i,9) = { 100e3*psia };
   wconinj(~i, 9) = cellfun(@(x) convertFrom(x, u.press), ...
                            wconinj(~i,9), 'UniformOutput', false);

   % General case
   c = [6                , 7, 10     ];
   u = [u.liqvol_r/u.time, 1, u.press];

   for n = 1 : numel(c),
      wconinj(:, c(n)) = cellfun(@(x) convertFrom(x, u(n)), ...
                                 wconinj(:, c(n)), 'UniformOutput', false);
   end
end

%--------------------------------------------------------------------------

function wconinj = convertWConInje(wconinj, u)
   % Surface flow rate
   is_gas = strcmpi(wconinj(:, 2), 'gas');
   is_oil = strcmpi(wconinj(:, 2), 'oil');

   r_unit         = repmat(u.liqvol_s / u.time, [numel(is_gas), 1]);
   r_unit(is_gas) = u.gasvol_s / u.time;

   wconinj(:, 5)  = arrayfun(@convertFrom,                   ...
                             vertcat(wconinj{:, 5}), r_unit, ...
                             'UniformOutput', false);

   % Vapourised gas concentration
   if any(is_gas | is_oil),
      R_unit = repmat(u.liqvol_s / u.gasvol_s, [numel(is_gas), 1]);
      R_unit(is_oil) = u.gasvol_s / u.liqvol_s;

      wconinj(:, 10) = arrayfun(@convertFrom,                    ...
                                vertcat(wconinj{:, 10}), R_unit, ...
                                'UniformOutput', false);
   end

   % General case
   c = [6                , 7      , 8      ];
   u = [u.liqvol_r/u.time, u.press, u.press];

   for n = 1 : numel(c),
      wconinj(:, c(n)) = cellfun(@(x) convertFrom(x, u(n)), ...
                                 wconinj(:, c(n)), 'UniformOutput', false);
   end
end

%--------------------------------------------------------------------------


function wconinj = convertWConInjh(wconinj, u)
   % Surface flow rate
   is_gas = strcmpi(wconinj(:, 2), 'gas');
   is_oil = strcmpi(wconinj(:, 2), 'oil');

   r_unit         = repmat(u.liqvol_s / u.time, [numel(is_gas), 1]);
   r_unit(is_gas) = u.gasvol_s / u.time;

   wconinj(:, 4)  = arrayfun(@convertFrom,                   ...
                             vertcat(wconinj{:, 4}), r_unit, ...
                             'UniformOutput', false);

   % Vapourised gas concentration
   if any(is_gas | is_oil),
      R_unit = repmat(u.liqvol_s / u.gasvol_s, [numel(is_gas), 1]);
      R_unit(is_oil) = u.gasvol_s / u.liqvol_s;

      wconinj(:, 8) = arrayfun(@convertFrom,                    ...
                                vertcat(wconinj{:, 8}), R_unit, ...
                                'UniformOutput', false);
   end

   % General case
   c = [5      , 6      ];
   u = [u.press, u.press];

   for n = 1 : numel(c),
      wconinj(:, c(n)) = cellfun(@(x) convertFrom(x, u(n)), ...
                                 wconinj(:, c(n)), 'UniformOutput', false);
   end
end

%--------------------------------------------------------------------------

function wcp = convertWConProd(wcp, u)
   resv = u.liqvol_r / u.time;
   lrat = u.liqvol_s / u.time;
   grat = u.gasvol_s / u.time;

   c = [4   , 5   , 6   , 7   , 8   , 9      , 10     ];
   u = [lrat, lrat, grat, lrat, resv, u.press, u.press];

   for n = 1 : numel(c),
      wcp(:, c(n)) = cellfun(@(x) convertFrom(x, u(n)), ...
                             wcp(:, c(n)), 'UniformOutput', false);
   end
end

%--------------------------------------------------------------------------

function wcp = convertWPolymer(wcp, u)

   c = [2        , 3        ];
   u = [u.density, u.density];

   for n = 1 : numel(c),
      wcp(:, c(n)) = cellfun(@(x) convertFrom(x, u(n)), ...
                             wcp(:, c(n)), 'UniformOutput', false);
   end
end
%--------------------------------------------------------------------------

function wcp = convertGConProd(wcp, u)
   resv = u.liqvol_r / u.time;
   lrat = u.liqvol_s / u.time;
   grat = u.gasvol_s / u.time;

   % Skip "calorific rate target or upper limit" (item 17).
   % We're isothermal...
   c = [3   , 4   , 5   , 6   , 14  , 15  , 16  , 17 , 18  , 19   ];
   u = [lrat, lrat, grat, lrat, resv, resv, grat, inf, grat, lrat ];

   for n = 1 : numel(c),
      wcp(:, c(n)) = cellfun(@(x) convertFrom(x, u(n)), ...
                             wcp(:, c(n)), 'UniformOutput', false);
   end
end

%--------------------------------------------------------------------------

function gconinj = convertGconInje(gconinj, u)
   % Surface flow rate
   is_gas = strcmpi(gconinj(:, 2), 'gas');

   r_unit         = repmat(u.liqvol_s / u.time, [numel(is_gas), 1]);
   r_unit(is_gas) = u.gasvol_s / u.time;

   gconinj(:, 4)  = arrayfun(@convertFrom,                   ...
                             vertcat(gconinj{:, 4}), r_unit, ...
                             'UniformOutput', false);

   % General case
   c = [5         , 13        ];
   u = [u.liqvol_r, u.gasvol_s] ./ u.time;

   for n = 1 : numel(c),
      gconinj(:, c(n)) = cellfun(@(x) convertFrom(x, u(n)), ...
                                 gconinj(:, c(n)), 'UniformOutput', false);
   end
end

%--------------------------------------------------------------------------

function wspecs = convertWellSpecs(wspecs, u)
   c = [5       , 7       ];
   u = [u.length, u.length];

   for n = 1 : numel(c),
      wspecs(:, c(n)) = cellfun(@(x) convertFrom(x, u(n)), ...
                                wspecs(:, c(n)), 'UniformOutput', false);
   end
end

%--------------------------------------------------------------------------

function s = id(s)
   s = ['convertDeckUnits:', s];
end
