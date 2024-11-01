function deck = readSOLUTION(fid, dirname, deck)

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

   [nc, ntequil, naqcon] = get_dimensions(deck);

   [sln, miss_kw] = get_state(deck);

   to_double = @(s) sscanf(regexprep(s, '[dD]', 'e'), '%f');

   kw = getEclipseKeyword(fid);
   in_section = ischar(kw);
   while in_section,
      switch kw,
         case 'AQUCT',
            tmpl(1 : 13) = { 'NaN' };
            tmpl(  9   ) = { '360' };

            data = readDefaultedKW(fid, tmpl, 'NRec', naqcon);

            sln.(kw) = cellfun(to_double, data);            clear data tmpl

         case 'AQUANCON',
            tmpl = [ repmat({ 'NaN' }, [1, 7]), ...
                     { '', '1', '-1', 'NO' } ];         % Negative => unset

            cols = [1:7, 9:10];
            data = readDefaultedKW(fid, tmpl, 'NRec', naqcon);

            data(:, cols) = cellfun(to_double, data(:, cols), ...
                                    'UniformOutput', false);

            sln.(kw) = data;                           clear data cols tmpl

         case 'BOX',
            boxKeyword(fid);

         case 'ENDBOX',
            endboxKeyword;

         case 'EQUIL',
            tmpl     = cell([1, 11]);
            tmpl(:)  = { '0' };
            tmpl{10} =   '1'  ;

            equil = readDefaultedKW(fid, tmpl, 'NRec', ntequil);
            sln.(kw) = cellfun(to_double, equil);                clear tmpl

         case 'DATUM',
            s = readRecordString(fid);
            data = to_double(s);

            if numel(data) ~= 1,
               error(['DATUM keyword must be followed by a single '  , ...
                      'real number (the datum depth for calculation ', ...
                      'of depth corrected pressures).']);
            end

            sln.(kw) = data;

         case { 'PBVD', 'PDVD', 'RSVD', 'RVVD' },
            for reg = 1 : ntequil,
               s             = readRecordString(fid);
               sln.(kw){reg} = reshape(to_double(s), 2, []) .';
            end

         case {'PRESSURE', 'RS', 'RV', 'SGAS', 'SOIL', 'SWAT'},
            sln.(kw) = readVector(fid, kw, nc);

         case {'ADD', 'COPY', 'EQUALS', 'MAXVALUE', ...
               'MINVALUE', 'MULTIPLY'},
            sln = applyOperator(sln, fid, kw);

         case {'ECHO', 'NOECHO'},
            kw = getEclipseKeyword(fid);
            continue;  % Ignore.  Not handled in MRST

         %-----------------------------------------------------------------
         % Sectioning keywords below.  Modifies flow of control.
         % Don't change unless absolutely required...
         %
         case 'END',
            % Logical end of input deck.
            %
            in_section   = false;
            deck.SOLUTION = sln;

         case 'SUMMARY',
            % Read next section (i.e., 'SOLUTION' -> 'SUMMARY' -> 'SCHEDULE')
            in_section = false;

            deck = set_state(deck, sln, miss_kw);

            deck = readSUMMARY(fid, dirname, deck);

         case 'SCHEDULE',
            % Read next section (i.e., 'SOLUTION' -> 'SCHEDULE', no 'SUMMARY')
            in_section = false;

            deck = set_state(deck, sln, miss_kw);

            % Restore default input box at end of section
            gridBox(defaultBox);

            deck = readSCHEDULE(fid, dirname, deck);

         case 'INCLUDE',
            % Handle 'INCLUDE' (recursion).
            deck = set_state(deck, sln, miss_kw);

            deck = readEclipseIncludeFile(@readSOLUTION, fid, dirname, deck);

            % Prepare for additional reading.
            [sln, miss_kw] = get_state(deck);

         otherwise,
            if ischar(kw),
               miss_kw = [ miss_kw, { kw } ];  %#ok
            end
      end

      % Get next keyword.
      kw = getEclipseKeyword(fid);
      in_section = in_section && ischar(kw);
   end

   deck = set_state(deck, sln, miss_kw);
end

%--------------------------------------------------------------------------

function [nc, ntequil, naqcon] = get_dimensions(deck)
   assert (isstruct(deck) && isfield(deck, 'RUNSPEC') && ...
           isstruct(deck.RUNSPEC));

   dims = reshape(deck.RUNSPEC.cartDims, 1, []);
   nc   = prod(dims);                 % Number of cells

   ntequil = 1;
   if isfield(deck.RUNSPEC, 'EQLDIMS'),
      ntequil = deck.RUNSPEC.EQLDIMS(1);
   end

   naqcon = 1;
   if isfield(deck.RUNSPEC, 'AQUDIMS'),
      naqcon = deck.RUNSPEC.AQUDIMS(2);
   end
   assert (ntequil >= 1);
end

%--------------------------------------------------------------------------

function [sln, miss_kw] = get_state(deck)
   sln     = deck.SOLUTION;
   miss_kw = deck.UnhandledKeywords.SOLUTION;
end

%--------------------------------------------------------------------------

function deck = set_state(deck, sln, miss_kw)
   deck.SOLUTION                   = sln;
   deck.UnhandledKeywords.SOLUTION = unique(miss_kw);
end