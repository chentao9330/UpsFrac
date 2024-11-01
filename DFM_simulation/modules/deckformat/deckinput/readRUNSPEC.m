function deck = readRUNSPEC(fid, dirname, deck)

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

   [rspec, miss_kw] = get_state(deck);

   double_conv = @(c) sscanf(regexprep(c, '[dD]', 'e'), '%f');
   to_double = @(data) cellfun(double_conv, data);

   kw = getEclipseKeyword(fid);
   in_section = ischar(kw);
   while in_section,
      switch kw,
         case 'AQUDIMS',
            tmpl = { '1', '1', '1', '36', '1', '1', '0', '0' };
            data = readDefaultedRecord(fid, tmpl);
            rspec.(kw) = to_double(data);  clear tmpl data

         case 'DIMENS',
            s = readRecordString(fid);
            rspec.cartDims = reshape(sscanf(s, '%f', 3), 1, []);
            rspec.DIMENS   = rspec.cartDims;

            % Set default input box corresponding to entire model.
            defaultBox(rspec.DIMENS);

         case 'ENDSCALE',
            tmpl        = { 'NODIR', 'REVERS', '1', '20', '0' };
            data        = readDefaultedRecord(fid, tmpl);
            data(3:end) = num2cell(to_double(data(3:end)));   clear tmpl
            rspec.(kw)  = data;

         case 'EQLDIMS',
            tmpl = {'1', '100', '50', '1', '50'};
            data = readDefaultedRecord(fid, tmpl);
            rspec.(kw) = to_double(data);  clear tmpl

         case 'GRIDOPTS',
            tmpl = { 'NO', '0', '0' };
            data = readDefaultedRecord(fid, tmpl);
            data(2:end) = cellfun(double_conv, data(2:end), ...
                                  'UniformOutput', false);
            rspec.(kw) = data;  clear tmpl data

         case 'MISCIBLE',
            tmpl = { '1', '20', 'NONE' };
            data = readDefaultedRecord(fid, tmpl);
            data(1:2) = cellfun(double_conv, data(1:2), ...
                                'UniformOutput', false);
            rspec.(kw) = data;  clear tmpl

         case 'REGDIMS',
            tmpl = { '1', '1', '0', '0', '0', '1', '0', '0', '0' };
            data = readDefaultedRecord(fid, tmpl);
            rspec.(kw) = to_double(data);  clear data tmpl

         case 'ROCKCOMP',
            tmpl = { 'REVERS', '1', 'NO', '', '0' };
            data = readDefaultedRecord(fid, tmpl);
            data([2,end]) = cellfun(double_conv, data([2,end]), ...
                                    'UniformOutput', false);
            rspec.(kw) = data;  clear tmpl data

         case 'START',
            s = readRecordString(fid);
            s = strrep(s, 'JLY', 'JUL');
            rspec.START = datenum(s, 'dd mmm yyyy');

         case 'SMRYDIMS',
            rspec.(kw) = readVector(fid, kw, 1);  % Or just ignored...

         case 'TABDIMS',
            tmpl = {'1', '1', '20', '20',  '1', '20', '20', '1', ...
                    '1', '1', '10',  '1', '-1',  '0',  '1'};
            data = readDefaultedRecord(fid, tmpl);
            rspec.(kw) = to_double(data);  clear tmpl

         case 'TITLE',
            rspec.(kw) = strtrim(fgetl(fid));

         case 'UDADIMS',
            tmpl = { '0', '0', '100' };
            data = readDefaultedRecord(fid, tmpl);
            rspec.(kw) = to_double(data);  clear tmpl data

         case 'UDQDIMS',
            tmpl = [{ '16', '16' }, repmat({ '0' }, [1, 8]), { 'N' }];
            data = readDefaultedRecord(fid, tmpl);
            data(1:end-1) = cellfun(double_conv, data(1:end-1), ...
                                       'UniformOutput', false);
            rspec.(kw) = data;   clear tmpl data

         case 'VFPPDIMS',
            tmpl(1:6) = { '1' };
            data = readDefaultedRecord(fid, tmpl);
            rspec.(kw) = to_double(data);  clear tmpl data

         case 'WELLDIMS',
            tmpl = {'10', 'inf', '5', '10', '5', '10', ...
                     '5',   '4', '3',  '0', '1',  '1'};
            data = readDefaultedRecord(fid, tmpl);  clear tmpl

            rspec.(kw) = to_double(data);
            if ~isfinite(rspec.(kw)(2)),
               rspec.(kw)(2) = rspec.cartDims(3);
            end

         case {'NOGRAV', 'IMPES',        ...
               'METRIC', 'FIELD', 'LAB', ...
               'WATER' , 'OIL'  , 'GAS', ...
               'DISGAS', 'VAPOIL',       ...
               'POLYMER', 'BRINE'},
            rspec.(kw) = true;

         case {'ECHO', 'NOECHO'},
            kw = getEclipseKeyword(fid);
            continue;  % Ignore.  Not handled in MRST

         %-----------------------------------------------------------------
         % Sectioning keywords below.  Modifies flow of control.
         % Don't change unless absolutely required...
         %
         case 'END',
            % Logical end of input deck.
            % Quite unusual (but nevertheless legal) in RUNSPEC.
            %
            in_section   = false;
            deck.RUNSPEC = rspec;

            % Restore default input box at end of section
            gridBox(defaultBox);

         case 'GRID',
            % Read next section (always 'GRID').
            in_section = false;

            deck = set_state(deck, rspec, miss_kw);

            % Restore default input box at end of section
            gridBox(defaultBox);

            deck = readGRID(fid, dirname, deck);

         case 'INCLUDE',
            % Handle 'INCLUDE' (recursion).
            deck = set_state(deck, rspec, miss_kw);

            deck = readEclipseIncludeFile(@readRUNSPEC, fid, dirname, deck);

            % Prepare for additional reading.
            [rspec, miss_kw] = get_state(deck);

         otherwise,
            if ischar(kw),
               miss_kw = [ miss_kw, { kw } ];  %#ok
            end
      end

      % Get next keyword.
      kw = getEclipseKeyword(fid);
      in_section = in_section && ischar(kw);
   end

   deck = set_state(deck, rspec, miss_kw);
end

%--------------------------------------------------------------------------

function [rspec, miss_kw] = get_state(deck)
   rspec   = deck.RUNSPEC;
   miss_kw = deck.UnhandledKeywords.RUNSPEC;
end

%--------------------------------------------------------------------------

function deck = set_state(deck, rspec, miss_kw)
   deck.RUNSPEC                   = rspec;
   deck.UnhandledKeywords.RUNSPEC = unique(miss_kw);
end
