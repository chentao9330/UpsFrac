function [press, s, rs, rv] = direct_assignment(G, deck, fluid)
%Initialisation through direct assignment of initial conditions.
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

   nm  = {'WATER', 'OIL' , 'GAS' };
   sat = {'SWAT' , 'SOIL', 'SGAS'};
   a   = isfield(deck.RUNSPEC , nm );
   i   = isfield(deck.SOLUTION, sat);

   check_consistency(deck, fluid, nm, sat, a, i)

   s     = initial_saturation  (G, deck, sat, a, i);
   press = initial_pressure    (G, deck, fluid, s, a);
   rs    = initial_solution_gor(G, deck, press, a);
   rv    = initial_vapour_ogr  (G, deck, press, a);
end

%--------------------------------------------------------------------------

function check_consistency(~, fluid, nm, sat, a, i)

   assert (all(strcmpi(nm(a), fluid.names)), ...
          ['Unexpected inconsistency in fluid phase names ', ...
           'versus declared phases.']);

   assert (any(sum(i) == sum(a) - [0, 1]), ...
           'Initial saturations may omit at most one phase.');

%{
   ngc = prod(deck.GRID.cartDims);
   assert (all(cellfun(@(f) numel(deck.SOLUTION.(f)), sat(i)) == ngc), ...
           'Initial saturations must be specified in each global cell.');
%}

   inconsist = i & ~a;   % Initial saturations and non-active phases.
   if any(inconsist),
      nl  = sprintf('\n');

      msg = ['Deck identifies initial saturations', nl, ' -'];
      msg = [msg, sprintf(' %s', sat{inconsist}), nl];
      msg = [msg, ...
             'but fails to declare corresponding active phases', nl, ' -'];
      msg = [msg, sprintf(' %s', nm{inconsist})];

      error(msg);
   end
end

%--------------------------------------------------------------------------

function s = initial_saturation(G, deck, sat, a, i)
   p  = cumsum(double(a));

   if isfield(G.cells, 'indexMap'),
      s0 = cellfun(@(f) deck.SOLUTION.(f)(G.cells.indexMap), ...
                   sat(i), 'UniformOutput', false);
   else
      s0 = cellfun(@(f) deck.SOLUTION.(f), sat(i), 'UniformOutput', false);
   end

   s          = nan([G.cells.num, p(end)]);
   s(:, p(i)) = [ s0{:} ];

   i = isnan(s(1, :));
   if any(i),
      s(:,i) = 1 - sum(s(:,~i), 2);
   end

   assert (all(all(isfinite(s))), ...
           'Initial saturation is non-finite?');

   assert (~ any(any((s < 0) | (s > 1))), ...
           'Unphysical initial saturations outside [0,1].');
end

%--------------------------------------------------------------------------

function press = initial_pressure(G, deck, fluid, s, a)
   if isfield(G.cells, 'indexMap'),
      imap = G.cells.indexMap;
   else
      imap = (1 : G.cells.num) .';
   end

   if isfield(deck.SOLUTION, 'PRESSURE'),

      po = reshape(deck.SOLUTION.PRESSURE(imap), [], 1);

   elseif isfield(deck.SOLUTION, 'PRVD'),

      if isfield(deck, 'REGIONS') && isfield(deck.REGIONS, 'EQLNUM'),
         eqlnum = reshape(deck.REGIONS.EQLNUM(imap), [], 1);
      else
         eqlnum = ones([G.cells.num, 1]);
      end

      po = nan([G.cells.num, 1]);

      for eqreg = reshape(find(accumarray(eqlnum, 1) > 0), 1, []),
         PRVD  = deck.SOLUTION.PRVD{eqreg};
         cells = eqlnum == eqreg;

         po(cells) = interp1(PRVD(:,1), PRVD(:,2), ...
                             G.cells.centroids(cells, 3));
      end

   else

      error(['Input deck does not define a means of obtaining ', ...
             'initial pressure distribution.']);

   end

   assert (all(isfinite(po)), ...
           'Some cells have non-finite initial pressures.');

   kr    = fluid.relperm(s);
   p     = cumsum(double(a));
   press = nan([G.cells.num, 1]);

   if a(2),                  % Active oil phase.
      i = kr(:, p(2)) > 0;   % Mobile oil phase.  Meaningful oil pressure.

      press(i) = po(i);
   end

   i = isnan(press);
   if any(i),
      [r, c] = find(kr(i,:) > 0);

      if isfield(fluid, 'pc'),
         pc = fluid.pc(s(i,:));
      else
         pc = zeros([numel(r), size(kr,2)]);
      end

      if a(2),
         pos = cumsum(double(a));
         assert (~ any(abs(pc(:, pos(2))) > 0), ...
                 'Non-zero oil capillary pressure?');
      end

      % Pcow = Po - Pw -> Pw = Po - Pcow
      press(i) = po(i) - pc(accumarray(r, c, [], @min));
   end
end

%--------------------------------------------------------------------------

function rs = initial_solution_gor(G, deck, press, a)

   if ~(a(2) && a(3)) || ~isfield(deck.PROPS, 'PVTO'),
      % Run that does not declare both OIL *and* GAS active (or a live oil
      % phase).  Set RS=0 and hope for the best.  Give a warning, though.

      rs = zeros([G.cells.num, 1]);

      if mrstVerbose,
         warning('ZeroRS:Unconditional', ...
                 'Assigning zero initial solution GOR (aka ''RS'').');
      end

   elseif isfield(deck.SOLUTION, 'RS'),

      rs = reshape(deck.SOLUTION.RS(G.cells.indexMap), [], 1);

   elseif isfield(deck.SOLUTION, 'PBUB'),

      [pvtnum, regno] = get_pvtnum(deck);

      rs   = nan([G.cells.num, 1]);
      pbub = reshape(deck.SOLUTION.PBUB(G.cells.indexMap), [], 1);

      for reg = regno,
         irs = interp_rs(deck.PROPS.PVTO{reg});

         cells     = pvtnum == reg;
         rs(cells) = irs(min(pbub(cells), press(cells)));
      end

      assert (all(isfinite(rs)), ...
              'Some cells not covered by Solution GOR initialisation.');

   else

      error('RS:InitUnknown', ...
            'Initial Solution GOR specified by unknown means.');

   end
end

%--------------------------------------------------------------------------

function rv = initial_vapour_ogr(G, deck, press, a)

   if ~(a(2) && a(3)) || ~isfield(deck.PROPS, 'PVTG'),
      % Run that does not declare both OIL *and* GAS active (or a wet gas
      % phase).  Set RV=0 and hope for the best.  Give a warning, though.

      rv = zeros([G.cells.num, 1]);

      if mrstVerbose,
         warning('ZeroRS:Unconditional', ...
                 'Assigning zero initial solution GOR (aka ''RS'').');
      end

   elseif isfield(deck.SOLUTION, 'RV'),

      rv = reshape(deck.SOLUTION.RV(G.cells.indexMap), [], 1);

   elseif isfield(deck.SOLUTION, 'PDEW'),

      [pvtnum, regno] = get_pvtnum(deck);

      rv     = nan([G.cells.num, 1]);
      pdew   = reshape(deck.SOLUTION.PDEW(G.cells.indexMap), [], 1);

      for reg = regno,
         irv = interp_rv(deck.PVTG{reg});

         cells     = pvtnum == reg;
         rv(cells) = irv(min(pdew(cells), press(cells)));
      end

   else

      error('RV:InitUnknown', ...
            'Initial Vapour OGR specified by unknown means.');

   end
end

%--------------------------------------------------------------------------

function [pvtnum, regno] = get_pvtnum(deck)
   if isfield(deck, 'REGIONS') && isfield(deck.REGIONS, 'PVTNUM'),
      pvtnum = reshape(deck.REGIONS.PVTNUM(G.cells.indexMap), [], 1);
   else
      pvtnum = ones([G.cells.num, 1]);
   end

   regno = reshape(find(accumarray(pvtnum, 1) > 0), 1, []);
end

%--------------------------------------------------------------------------

function irs = interp_rs(PVTO)
   irs = @(p) interp1(PVTO.data(PVTO.pos(1:end-1), 1), ...
                      PVTO.key                       , ...
                      p, 'linear', 'extrap');
end

%--------------------------------------------------------------------------

function irv = interp_rv(PVTG)
   irv = @(p) interp1(PVTG.key                       , ...
                      PVTG.data(PVTG.pos(1:end-1), 1), ...
                      p, 'linear', 'extrap');

end
