function require(varargin)
%Announce and enforce module dependency.
%
% SYNOPSIS:
%   require modulename
%   require('modulename')
%
% PARAMETERS:
%   modulename - string.
%
% SEE ALSO:
%   mrstModule

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

   modules = @(lst) cellfun(@fileparts, lst, 'UniformOutput', false);

   modlist          = mrstModule;
   [ignore, active] = modules(modlist);   %#ok
   [ignore, needed] = modules(varargin);  %#ok

   if ~isempty(active),
      % 'Missing' is the subset of 'needed' modules @NOT in 'active'.
      missing = subset_op(needed, active, @not);
   else
      missing = needed;
   end

   if ~isempty(missing),
      [ignore, i, ignore] = unique(missing);  %#ok
      missing = missing(sort(i)); % SORT by *LAST* occurrence of ind. module

      plural = '';
      if numel(missing) ~= 1, plural = 's'; end

      e = sprintf('Did you forget to include the %s module%s?', ...
                  quote_list(missing), plural);

      % 'm' includes '.', and '..', but that's not a problem in this
      % particular case.
      %
      m = dir(fullfile(ROOTDIR, 'modules'));

      % 'Known' is the subset of 'missing' modules in module directory.
      known   = subset_op(missing, { m([m.isdir]).name }, @(x) x);
      plural  = ' is';

      if numel(known) > 1, plural = 's are'; end

      if ~isempty(known),
         a = sprintf(['The %s module%s known to MRST and may be ', ...
                      'activated using\n\n\tmrstModule add%s'], ...
                      quote_list(known), plural, stringify_list(known));

         if numel(known) < numel(missing),
            plural = '';
            if numel(known) < numel(missing) - 1, plural = 's'; end

            a = [a, sprintf('\n\nThe remaining module%s', plural), ...
                 ' ', hercule, '.'];
         end
      else
         plural1 = '';
         plural2 = ' is';
         if numel(missing) > 1,
            plural1 = 'se';
            plural2 = 's are';
         end

         a = [sprintf('The%s module%s not known to MRST and ', ...
              plural1, plural2),  hercule, sprintf('\n'), ...
              'before activation through ''mrstModule add''.'];
      end

      throwAsCaller(MException('MRST:MissingModule', '%s\n%s',e, a));
   end
end

%--------------------------------------------------------------------------

function v = subset_op(needed, known, op)
% Slightly generalised SETDIFF; preserves original ordering of 'needed'.

   nn = numel(needed);

   [i, j] = ndgrid(1:nn, 1:numel(known));
   t      = strcmp(reshape(needed(i(:)), [], 1), ...
                   reshape(known (j(:)), [], 1));

   present       = false([nn, 1]);
   present(i(t)) = true;

   v = needed(op(present));
end

%--------------------------------------------------------------------------

function s = quote_list(v)
   assert (iscellstr(v), 'List must be cell array of strings');

   if numel(v) == 1,
      s = ['''', v{1}, ''''];
   elseif numel(v) > 1,
      if numel(v) > 2,
         sep = ', ';
      else
         sep = ' ';
      end

      s = [sprintf(['''%s''', sep], v{1:end-1}), ...
           'and ''', v{end}, ''''];
   else
      s = '';
   end
end

%--------------------------------------------------------------------------

function s = stringify_list(v)
   assert (iscellstr(v), 'List must be cell array of strings');

   if numel(v) == 1,
      s = [' ', v{1}];
   elseif numel(v) > 1,
      s = sprintf(' %s', v{:});
   else
      s = '';
   end
end

%--------------------------------------------------------------------------

function s = hercule
   s = 'must be searched for by method and reason';
end
