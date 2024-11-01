function [press, s, rs, rv] = equilibration(G, deck, fluid)
%Equilibration facility to handle EQUIL keyword.
%
% Internal helper function.  Interface subject to change.
%
% SEE ALSO:
%   initEclipseState.

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

   assert (isfield(G.cells, 'centroids'), ...
           'Input grid must be equipped with geometric information.');

   if isfield(deck, 'REGIONS') && isfield(deck.REGIONS, 'EQLNUM'),
      eqlnum = reshape(deck.REGIONS.EQLNUM(G.cells.indexMap), [], 1);
   else
      eqlnum = ones([G.cells.num, 1]);
   end

   if isfield(deck, 'REGIONS') && isfield(deck.REGIONS, 'PVTNUM'),
      pvtnum = reshape(deck.REGIONS.PVTNUM(G.cells.indexMap), [], 1);
   else
      pvtnum = ones([G.cells.num, 1]);
   end

   np        = numel(fluid.names);
   pressure  = zeros([G.cells.num, np]);
   rs        = zeros([G.cells.num, 1 ]);
   rv        = zeros([G.cells.num, 1 ]);

   for eqreg = reshape(find(accumarray(eqlnum, 1) > 0), 1, []),
      cells  = eqlnum == eqreg;
      pvtreg = find(accumarray(pvtnum(cells), 1) > 0);

      if numel(pvtreg) > 1,
         error(id('EquilReg:InconsistentPVTReg')            , ...
              ['Cells in equilibration region %d reference ', ...
               '%d PVT tables.'], eqreg, numel(pvtreg));
      end

      [pressure(cells, :), rs(cells, :), rv(cells, :)] = ...
         equilibrate_region(deck, G, cells, fluid, eqreg, pvtreg);
   end

   s  = initsat(deck, pressure, G.cells.centroids(:,3), G.cells.indexMap);
   kr = fluid.relperm(s);

   % Don't pretend there's a meaningful phase pressure unless the
   % corresponding phase is mobile.
   pressure(~ (kr > 0.0)) = nan;

   liq      = strcmpi(fluid.names, 'oil');
   i        = isfinite(pressure(:, liq));
   press    = zeros([G.cells.num, 1]);
   press(i) = pressure(i, liq);
end

%--------------------------------------------------------------------------
% Private helpers follow.
%--------------------------------------------------------------------------

function [press, rs, rv] = ...
      equilibrate_region(deck, G, cells, fluid, eqreg, pvtreg)
   equil = deck.SOLUTION.EQUIL(eqreg, :);

   Z0   = equil(1);   P0 = equil(2);
   Zwoc = equil(3);
   Zgoc = equil(5);

   Zc   = G.cells.centroids(cells, 3);
   Zmin = min(Zc);  Zmax = max(Zc);

   Rs = define_dissolution(deck, Z0, P0, Zgoc, eqreg, pvtreg);
   Rv = define_evaporation(deck, Z0, P0, Zgoc, eqreg, pvtreg);

   if ~((Zgoc > Z0) || (Z0 > Zwoc))
      % Datum depth in oil zone  (Zgoc <= Z0 <= Zwoc)
      press = equilibrate_OWG(Zc, [Zmin, Zmax], fluid, Rs, Rv, equil(1:6));
   else
      error('Datum not in oil zone')
   end

   rs = Rs(Zc, press(:,2)); % min(Rs(z), Rs(p_oil))
   rv = Rv(Zc, press(:,3)); % min(Rv(z), Rv(p_gas))
end

%--------------------------------------------------------------------------

function s = initsat(deck, press, Zc, indexMap)
   % Invert capillary pressure function to derive initial saturations.

   if isfield(deck.REGIONS, 'SATNUM'),
      satnum = deck.REGIONS.SATNUM(indexMap);
   else
      satnum = ones([numel(indexMap), 1]);
   end

   Pcow = press(:,2) - press(:,1);   % P_o - P_w
   Pcog = press(:,3) - press(:,2);   % P_g - P_o

   s = zeros(size(press));

   [inv_Pcow, inv_Pcog] = pc_inv(deck);

   for reg = reshape(find(accumarray(satnum, 1) > 0), 1, []),
      i = satnum == reg;

      s(i,1) = inv_Pcow{reg}(Pcow(i), Zc(i));
      s(i,3) = inv_Pcog{reg}(Pcog(i), Zc(i));
   end

   % Define oil saturation to fill pore space.
   s(:,2) = 1 - sum(s, 2);

   assert (~any(any(s < 0)), 'Negative saturations during equilibrium');
end

%--------------------------------------------------------------------------

function Rs = define_dissolution(deck, Z0, P0, Zgoc, eqreg, pvtreg)
   if isfield(deck.PROPS, 'PVTO'),
      PVTO = deck.PROPS.PVTO{pvtreg};

      interp_Rs = @(p) interp1(PVTO.data(PVTO.pos(1:end-1), 1), ...
                               PVTO.key                       , ...
                               p, 'linear', 'extrap');

      if     isfield(deck.SOLUTION, 'RSVD'),

         RSVD = deck.SOLUTION.RSVD{eqreg};
         Rs   = @(z, p) interp1(RSVD(:,1), RSVD(:,2), ...
                                z, 'linear', 'extrap');

      elseif isfield(deck.SOLUTION, 'PBVD'),

         PBVD = deck.SOLUTION.PBVD{eqreg};
         Rs   = @(z, p) interp_Rs(min(p, interp1(PBVD(:,1), PBVD(:,2), ...
                                                 z, 'linear', 'extrap')));

      else
         if Z0 ~= Zgoc,
            error(id('DefaultRS:DepthMismatch')                  , ...
                 ['Table %d: Datum depth must coincide with GOC ', ...
                  'in absence of explicit solubility data for '  , ...
                  'equilibration.'], eqreg);
         end

         Rsmax = interp_Rs(P0);
         Rs    = @(z, p) min(interp_Rs(p), Rsmax);
      end

   else
      % Immiscible
      Rs = @(varargin) 0;
   end
end

%--------------------------------------------------------------------------

function Rv = define_evaporation(deck, Z0, P0, Zgoc, eqreg, pvtreg)
   if isfield(deck.PROPS, 'PVTG'),
      PVTG = deck.PROPS.PVTG{pvtreg};

      interp_Rv = @(p) interp1(PVTG.key                       , ...
                               PVTG.data(PVTG.pos(1:end-1), 1), ...
                               p, 'linear', 'extrap');

      if     isfield(deck.SOLUTION, 'RVVD'),

         RVVD = deck.SOLUTION.RVVD{eqreg};
         Rv   = @(z, p) interp1(RVVD(:,1), RVVD(:,2), ...
                                z, 'linear', 'extrap');

      elseif isfield(deck.SOLUTION, 'PDVD'),

         PDVD = deck.SOLUTION.PDVD{eqreg};
         Rv   = @(z, p) interp_Rv(min(p, interp1(PDVD(:,1), PDVD(:,2), ...
                                                 z, 'linear', 'extrap')));

      else
         if Z0 ~= Zgoc,
            error(id('DefaultRV:DepthMismatch')                  , ...
                 ['Table %d: Datum depth must coincide with GOC ', ...
                  'in absence of explicit solubility data for '  , ...
                  'equilibration.'], eqreg);
         end

         Rv = @(z, p) interp_Rv(min(p, P0));
      end

   else
      % Immiscible
      Rv = @(varargin) zeros([numel(varargin{1}), 1]);
   end
end

%--------------------------------------------------------------------------

function press = equilibrate_OWG(Zc, bdry, fluid, Rs, Rv, equil)
   odeopts = odeset('AbsTol', 1.0e-10, 'RelTol', 5.0e-8);

   press   = zeros([numel(Zc), 3]);

   Z0   = equil(1);   P0       = equil(2);
   Zwoc = equil(3);   Pcow_woc = equil(4);  % Water-oil contact
   Zgoc = equil(5);   Pcog_goc = equil(6);  % Gas-oil contact

   assert (~(Zwoc < Zgoc), 'WOC above GOC?!?');

   Zmin = bdry(1);
   Zmax = bdry(2);

   O_up = [Z0  , Zmin];  O_down = [Z0  , Zmax];
   W_up = [Zwoc, Zmin];  W_down = [Zwoc, Zmax];
   G_up = [Zgoc, Zmin];  G_down = [Zgoc, Zmax];

   norm_g = norm(gravity());

   %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   % 1a) Oil pressure.
   %
   if ~(Z0 < Zmin) && abs(diff(O_up)) > 0,
      sol_u = ode45(@(z, p) dp_o(z, p, fluid, norm_g, Rs), ...
                    O_up, P0, odeopts);
   else
      sol_u = [];
   end

   if ~(Z0 > Zmax) && abs(diff(O_down)) > 0,
      sol_d = ode45(@(z, p) dp_o(z, p, fluid, norm_g, Rs), ...
                    O_down, P0, odeopts);
   else
      sol_d = [];
   end

   %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   % 1b) Derive phase pressures at WOC and GOC
   %
   if     diff(W_up  ) > 0,  % WOC above reservoir   (unexpected)
      P0w =  inf;
   elseif diff(W_down) < 0,  % WOC below reservoir
      P0w = -inf;
   else                      % WOC in    reservoir
      if Zwoc < Z0,
         assert (~isempty(sol_u), ...
                 'Internal error defining oil pressure equilibrium.');

         % Datum set in W-zone                       (unexpected)
         P0w = deval(sol_u, Zwoc) - Pcow_woc;
      else
         assert (~isempty(sol_d), ...
                 'Internal error defining oil pressure equilibrium.');

         % Datum set above W-zone
         P0w = deval(sol_d, Zwoc) - Pcow_woc;
      end
   end

   if     diff(G_up  ) > 0,  % GOC above reservoir
      P0g = -inf;
   elseif diff(G_down) < 0,  % GOC below reservoir   (unexpected)
      P0g =  inf;
   else                      % WOC in    reservoir
      if Zgoc < Z0,
         assert (~isempty(sol_u), ...
                 'Internal error defining oil pressure equilibrium.');

         % Datum set below G-zone
         P0g = deval(sol_u, Zgoc) - Pcog_goc;
      else
         assert (~isempty(sol_d), ...
                 'Internal error defining oil pressure equilibrium.');

         % Datum set in G-zone                       (unexpected)
         P0g = deval(sol_d, Zgoc) - Pcog_goc;
      end
   end

   cu = Zc < Z0;
   if any(cu),
      assert (~isempty(sol_u), 'Internal error.');
      press(cu, 2) = deval(sol_u, Zc( cu));
   end
   if ~all(cu),
      cdown = Zc > Z0;

      if any(cdown)
         assert (~isempty(sol_d), 'Internal error.');
         press(cdown, 2) = deval(sol_d, Zc(cdown));
      end

      press(~(cu | cdown), 2) = P0;
   end

   %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   % 2) Water pressure.
   %
   if isfinite(P0w),
      % Water present during initialisation/equilibration.

      sol_u = ode45(@(z, p) dp_w(p, fluid, norm_g), W_up  , P0w, odeopts);
      sol_d = ode45(@(z, p) dp_w(p, fluid, norm_g), W_down, P0w, odeopts);

      cu = Zc < Zwoc;
      if  any(cu), press( cu, 1) = deval(sol_u, Zc( cu)); end
      if ~all(cu), press(~cu, 1) = deval(sol_d, Zc(~cu)); end
   else
      press(:, 1) = P0w;
   end

   %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   % 3) Gas pressure.
   %
   if isfinite(P0g),
      % Gas present during initialisation/equilibration.

      sol_u = ode45(@(z, p) dp_g(z, p, fluid, norm_g, Rv), ...
                    G_up, P0g, odeopts);
      sol_d = ode45(@(z, p) dp_g(z, p, fluid, norm_g, Rv), ...
                    G_down, P0g, odeopts);

      cu = Zc < Zgoc;
      if  any(cu), press( cu, 3) = deval(sol_u, Zc( cu)); end
      if ~all(cu), press(~cu, 3) = deval(sol_d, Zc(~cu)); end
   else
      press(:, 3) = P0g;
   end
end

%--------------------------------------------------------------------------

function [iPcow, iPcog] = pc_inv(deck)
   if isfield(deck.PROPS, 'SWOF'),
      iPcow = cell([1, numel(deck.PROPS.SWOF)]);

      for reg = 1 : numel(iPcow),
         t    = deck.PROPS.SWOF{reg}(:, [1, end]);
         zwoc = deck.SOLUTION.EQUIL(reg, 3);

         iPcow{reg} = ...
            @(pcow, depth) inv_interp(pcow, depth, zwoc, t, true);
      end

   elseif isfield(deck.PROPS, 'SWFN'),
      iPcow = cell([1, numel(deck.PROPS.SWFN)]);

      for reg = 1 : numel(iPcow),
         t    = deck.PROPS.SWFN{reg}(:, [1, end]);
         zwoc = deck.SOLUTION.EQUIL(reg, 3);

         iPcow{reg} = ...
            @(pcow, depth) inv_interp(pcow, depth, zwoc, t, true);
      end

   else
      iPcow = cell([1, size(deck.SOLUTION.EQUIL, 1)]);

      for reg = 1 : numel(iPcow),
         assert (deck.SOLUTION.EQUIL(reg, 4) == 0, 'HuH?');
         zwoc = deck.SOLUTION.EQUIL(reg, 3);

         iPcow{reg} = @(pcow, depth) double(depth > zwoc);
      end

      warning('iPcow:Undefined', ...
             ['Inverse water capillary pressure curve undefined for ', ...
              'saturation functions other than ''SWOF'' or ''SWFN''.']);
   end

   if isfield(deck.PROPS, 'SGOF'),
      iPcog = cell([1, numel(deck.PROPS.SGOF)]);

      for reg = 1 : numel(iPcog),
         t    = deck.PROPS.SGOF{reg}(:, [1, end]);
         zgoc = deck.SOLUTION.EQUIL(reg, 5);

         iPcog{reg} = ...
            @(pcog, depth) inv_interp(pcog, depth, zgoc, t, false);
      end

   elseif isfield(deck.PROPS, 'SGFN'),
      iPcog = cell([1, numel(deck.PROPS.SGFN)]);

      for reg = 1 : numel(iPcog),
         t    = deck.PROPS.SGFN{reg}(:, [1, end]);
         zgoc = deck.SOLUTION.EQUIL(reg, 5);

         iPcog{reg} = ...
            @(pcog, depth) inv_interp(pcog, depth, zgoc, t, false);
      end

   else
      iPcog = cell([1, size(deck.SOLUTION.EQUIL, 1)]);

      for reg = 1 : numel(iPcog),
         assert (deck.SOLUTION.EQUIL(reg, 6) == 0, 'HuH?');
         zgoc = deck.SOLUTION.EQUIL(5);

         iPcog{reg}= @(pcow, depth) double(depth < zgoc);
      end

      warning('iPcog:Undefined', ...
             ['Inverse gas capillary pressure curve undefined for ', ...
              'saturation functions other than ''SGOF'' or ''SGFN''.']);
   end
end

%--------------------------------------------------------------------------

function s = id(s)
   s = ['EQUIL:', s];
end

%--------------------------------------------------------------------------

function dp = dp_o(z, p, fluid, norm_g, Rs)
   R = Rs(z, p);

   [rho, rho] = fluid.pvt(p, [0, 1, R]);  %#ok

   dp = rho(2) * norm_g;
end

%--------------------------------------------------------------------------

function dp = dp_w(p, fluid, norm_g)

   [rho, rho] = fluid.pvt(p, [1, 0, 0]);  %#ok

   dp = rho(1) * norm_g;
end

%--------------------------------------------------------------------------

function dp = dp_g(z, p, fluid, norm_g, Rv)
   R = Rv(z, p);

   [rho, rho] = fluid.pvt(p, [0, R, 1]);  %#ok

   dp = rho(3) * norm_g;
end

%--------------------------------------------------------------------------

function s = inv_interp(p, d, dc, tbl, is_increasing)
   assert (all(diff(tbl(:, 1)) > 0), ...
           'Saturation values are not monotonically increasing.');

   s = nan(size(p));

   if norm(diff(tbl(:,2))) > 0,
      [min_p, min_i] = min(tbl(:,2));
      [max_p, max_i] = max(tbl(:,2));

      i = ~((p < min_p) | (p > max_p));

      if max_p ~= min_p,
         s(i)      = interp1(tbl(:,2), tbl(:,1), p(i));
      end
      s(p < min_p) = tbl(min_i, 1);
      s(p > max_p) = tbl(max_i, 1);
   else
      % Constant capillary pressure function.  Use step definition.

      m = min(tbl(:,1));   % Minimum (connate) saturation
      M = max(tbl(:,1));   % Maximum saturation

      i = d < dc;

      if is_increasing,
         s( i) = m; s(~i) = M;
      else
         s(~i) = m; s( i) = M;
      end
   end

   assert (~any(isnan(s)), ...
           'Internal error during capillary pressure inversion.');
end
