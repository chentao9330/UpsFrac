function w  = readWellKW(fid, w, kw)
%Read well definitions from an ECLIPSE Deck specification.
%
% SYNOPSIS:
%   w = readWellKW(fid, w, kw)
%
% PARAMETERS:
%   fid - Valid file identifier as obtained from FOPEN.
%
%   w   - Structure array containing current well configuration.  May be
%         empty if the well keyword is 'WELSPECS'.
%
%   kw  - Current well specification keyword.
%
% RETURNS:
%   w   - Updated structure array containing new data for keyword 'kw'.
%
% COMMENTS:
%   The currently recognized well keywords are:
%      'WELSPECS', 'COMPDAT', 'WCONINJE', 'WCONPROD', and 'WCONHIST'.
%
% SEE ALSO:
%   readGRDECL, processWells.

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

   % Unless we're defining a new set of wells, there had better be some
   % previously defined wells.
   %
   assert (strcmp(kw, 'WELSPECS') || ~isempty(w.WELSPECS), ...
          ['Well keyword ''%s'' encountered before any wells have ', ...
           'been declared using ''WELSPECS''.'], kw);

   switch kw,
      case 'WELSPECS', w = readWellSpec(fid, w);
      case 'COMPDAT' , w = readCompDat (fid, w);
      case 'WCONHIST', w = readWConHist(fid, w);
      case 'WCONINJ' , w = readWConInj (fid, w);
      case 'WCONINJE', w = readWConInje(fid, w);
      case 'WCONINJH', w = readWConInjh(fid, w);
      case 'WCONPROD', w = readWConProd(fid, w);
      case 'GRUPTREE', w = readGrupTree(fid, w);
      case 'WGRUPCON', w = readWGrupCon(fid, w);
      case 'GCONPROD', w = readGConProd(fid, w);
      case 'WPOLYMER', w = readWPolymer(fid, w);
      otherwise
         fclose(fid);
         error(msgid('WellKW:Unsupported'), ...
               'Well keyword ''%s'' is not currently supported.', kw);
   end
end

%--------------------------------------------------------------------------
% Helpers follow.
%--------------------------------------------------------------------------

function w = readGrupTree(fid, w)
   template = {'Default', 'Default'};
   gruptree = readDefaultedKW(fid, template);
   w.GRUPTREE = [w.GRUPTREE; gruptree];
end

%--------------------------------------------------------------------------

function w = readWellSpec(fid, w)
   template = {'Default', 'Default', 'Default', 'Default', 'NaN', ...
               'Default', '0.0', 'STD', 'SHUT', 'NO', '0', 'SEG', '0'};
   numeric  = [3, 4, 5, 7, 11, 13];

   WellSpec = readDefaultedKW(fid, template);

   if ~isempty(WellSpec),
      WellSpec = toDouble(WellSpec, numeric);

      if ~isempty(w.WELSPECS),
         [tf, loc] = ismember(WellSpec(:,1), w.WELSPECS(:,1));
         ws = w.WELSPECS;
         if any(tf),
            ws(loc(tf),:) = WellSpec(tf,:);
         end
         WellSpec = [ws; WellSpec(~tf,:)];
      end
      w.WELSPECS = WellSpec;
   end
end

%--------------------------------------------------------------------------

function w = readCompDat(fid, w)
   template = {'Default', '-1', '-1', 'Default', 'Default', 'OPEN', ...
               '-1', '-1.0', '0.0', '-1.0', '0.0', 'Default', 'Z', '-1'};
   numeric  = [2:5, 7:11, 14];

   CompDat = readDefaultedKW(fid, template);

   % Deliberate misrepresentation of complexity.  Nevertheless a reasonable
   % approach in most (all?) of our cases...
   w.COMPDAT = [w.COMPDAT; toDouble(CompDat, numeric)];
end

%--------------------------------------------------------------------------

function w = readWConInje(fid, w)
   template = {'Default', 'Default', 'OPEN', 'Default', 'inf', 'inf', ...
               'NaN', 'inf', '0', '0.0'};
   numeric  = 5 : numel(template);

   data = readDefaultedKW(fid, template);
   data = toDouble(data, numeric);

   w.WCONINJE = appendSpec(w.WCONINJE, data, w.WELSPECS(:,1));
end

%--------------------------------------------------------------------------

function w = readWConInjh(fid, w)
   template = {'Default', 'Default', 'OPEN', '0.0', '0.0', '0.0', 'inf'...
               '0.0', '0', '0.0', '0.0'};
   numeric  = 4 : numel(template);

   data = readDefaultedKW(fid, template);
   data = toDouble(data, numeric);

   w.WCONINJH = appendSpec(w.WCONINJH, data, w.WELSPECS(:,1));
end

%--------------------------------------------------------------------------

function w = readWConInj(fid, w)
   template = {'Default', 'Default', 'OPEN', 'Default', 'inf', 'inf', ...
               '0.0', 'NONE', '-1', 'inf', '0.0', '0.0'};
   numeric  = [5:7, 9:numel(template)];

   data = readDefaultedKW(fid, template);
   data = toDouble(data, numeric);

   w.WCONINJ = appendSpec(w.WCONINJ, data, w.WELSPECS(:,1));
end

%--------------------------------------------------------------------------

function w = readWConProd(fid, w)
   template = {'Default', 'OPEN', 'Default', 'inf', 'inf', 'inf', ...
               'inf', 'inf', 'NaN', '0.0', '0', '0.0'};
   numeric  = 4 : numel(template);

   data = readDefaultedKW(fid, template);
   data = toDouble(data, numeric);

   w.WCONPROD = appendSpec(w.WCONPROD, data, w.WELSPECS(:,1));
end

%--------------------------------------------------------------------------

function w = readGConProd(fid, w)
   template = {'Default', 'NONE', 'Default', 'inf', 'inf', 'inf', ...
               'NONE', 'YES', 'inf', '', ...
               'NONE', ...
               'NONE', ...
               'NONE', 'inf', 'inf', 'inf', 'inf', 'inf', 'inf', ...
               'inf' , ...
               'NONE', ...
               };
   numeric  = [3:6, 9, 14:20];

   data = readDefaultedKW(fid, template);
   data = toDouble(data, numeric);

   w.GCONPROD = [w.GCONPROD; data];
end

%--------------------------------------------------------------------------

function w = readWGrupCon(fid, w)
   template = {'Default', 'YES', '-1.0', '', '1.0'};
   numeric  = [3,5];

   data = readDefaultedKW(fid, template);
   data = toDouble(data, numeric);

   w.WGRUPCON = [w.WGRUPCON; data];
end

%--------------------------------------------------------------------------

function w = readWConHist(fid, w)
   template = {'Default', 'OPEN', 'Default', '0.0', '0.0', '0.0', ...
               '0', 'Default', '0.0', '0.0', '0.0'};
   numeric  = [4 : 7, 9 : numel(template)];

   data = readDefaultedKW(fid, template);
   data = toDouble(data, numeric);

   w.WCONHIST = appendSpec(w.WCONHIST, data, w.WELSPECS(:,1));
end

%--------------------------------------------------------------------------

function w = readWPolymer(fid, w)
   template = {'Default', '0.0', '0.0', 'Default', 'Default'};
   numeric  = [2, 3];

   data = readDefaultedKW(fid, template);
   data = toDouble(data, numeric);

   w.WPOLYMER = appendSpec(w.WPOLYMER, data, w.WELSPECS(:,1));
end
%--------------------------------------------------------------------------

function T = toDouble(T, col)
   try
      T(:,col) = reshape(cellfun(@(v) sscanf(v, '%f'), T(:,col), ...
                                 'UniformOutput', false),        ...
                         [], numel(col));
   catch ME
      error('NonNumeric:Unexpected', ...
            'Unexpected non-numeric data encountered.');
   end
end

%--------------------------------------------------------------------------

function table = appendSpec(table, data, wells)
   matches = @(a,p) cellfun(@(c) ~isempty(c), regexp(a, p, 'match'));

   % Wildcard specs in new data?
   is_wc = matches(data(:,1), '\w+\*\s*$');
   for i = reshape(find(is_wc), 1, []),      % Row shape is essential here.
      patt = strrep(data{i,1}, '*', '.*');

      if ~isempty(table),
         % Remove any previous specs matching this 'patt'ern.
         table(matches(table(:,1), patt), :) = [];
      end

      match_wells = matches(wells, patt);
      nmatch      = sum(double(match_wells));
      if nmatch == 0,
         warning(msgid('Well:Unknown'), ...
                ['Well control wildcard ''%s'' does not match ', ...
                 'any previously defined wells.  Ignored.'], data{i,1});
      end
      append      = repmat(data(i,:), [nmatch, 1]);
      append(:,1) = wells(match_wells);

      table = [table; append];  %#ok       % Willfully ignore MLINT advice.
   end

   % Remove wildcard data.
   data(is_wc, :) = [];

   % Ensure that any remaining injection specs refer to known wells.
   if ~isempty(data) && ~all(ismember(data(:,1), wells)),
      error(msgid('Well:Unknown'), ...
           ['Well control specified in at least one undefined well.\n', ...
            'Attempted coup d''état foiled.']);
   end

   % Complete keyword spec for this configuration.
   table = [table; data];
end
