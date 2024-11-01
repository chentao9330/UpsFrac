function varargout = mrstModule(varargin)
%Query or modify list of activated add-on MRST modules
%
% SYNOPSIS:
%   Either of the modes
%      1) mrstModule <command> [module list]
%      2) modules = mrstModule
%
% PARAMETERS:
%   Mode 1)
%     <Command>     - One of the explicit verbs 'add', 'clear', 'list' or
%                     'reset'. The semantics of the command verbs are as
%                     follows:
%
%                       o) add   -- Activate specified modules from the
%                                   [module list].  Modules already
%                                   activated are moved to the beginning of
%                                   MATLAB's search path and remain active.
%
%                       o) clear -- Deactivate all modules.  An explicit
%                                   module list, if present, is ignored.
%
%                       o) list  -- Display list of currently active
%                                   modules in command window.  An explicit
%                                   module list, if present, is ignored.
%
%                       o) reset -- Convenience verb.  Equivalent to the
%                                   verb sequence:
%                                      mrstModule clear
%                                      mrstModule add [module list]
%
%     [module list] - A sequence of strings naming individual add-on
%                     modules for MRST.  A module string/name may be either
%                     of the following:
%
%                       o) A relative path-name such as 'agglomeration' or
%                          'eclipse/resultinput'.  The name is interpreted
%                          as a directory relative to the default MRST
%                          module directory,
%
%                              fullfile(ROOTDIR, 'modules')
%
%                          If no such directory exists, as determined by
%                          ISDIR, the requested module is ignored.
%
%                          The string may include the path-name component
%                          separator, FILESEP.
%
%                       o) An absolute path-name identifying an arbitrary
%                          directory on the local computer system.  If the
%                          directory does not exist, the requested module
%                          is ignored.
%
%   Mode 2)
%     None.
%
% RETURNS:
%   Mode 1)
%     Nothing.
%
%   Mode 2)
%     modules - List, represented as a cell array of strings, of the
%               currently active add-on modules.
%
% EXAMPLES:
%    mrstModule add eclipse
%    mrstModule add /exp/module
%    mrstModule add ../utility/module
%
%    mrstModule list
%    mrstModule clear
%
% SEE ALSO:
%   ROOTDIR, isdir, filesep.

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

   persistent MODLIST

   % Enforce cell array type.
   if isempty(MODLIST), MODLIST = {}; end

   if nargin > 0,
      assert (all(cellfun(@ischar, varargin)), ...
              'All parameters must be strings.');

      cmd  = varargin{1};
      mods = varargin(2 : end);

      switch lower(cmd),
         case 'add',
            MODLIST = prune_modules(MODLIST);

            if ~ (isempty(mods) || all(cellfun(@isempty, mods))),
               mlock
               MODLIST = add_modules(MODLIST, mods);
            else
               munlock
            end

         case 'clear',
            munlock
            MODLIST = clear_modules(MODLIST);

         case 'list',
            MODLIST = prune_modules(MODLIST);

            if numel(MODLIST) > 0,
               fprintf('Currently active MRST modules\n');
               fprintf('  * %s\n', MODLIST{:});
            else
               if mislocked, munlock, end;
               fprintf('No active MRST modules\n');
            end

         case 'reset',
            munlock
            MODLIST = clear_modules(MODLIST);

            if ~ (isempty(mods) || all(cellfun(@isempty, mods))),
               mlock
               MODLIST = add_modules(MODLIST, mods);
            end

         otherwise,
            error(msgid('Command:Unsupported'), ...
                 ['Command word ''%s'' unsupported. Must be one of ', ...
                  '''add'', ''clear'', ''list'', or ''reset''.'], cmd);
      end

   elseif nargout == 0,

      mrstModule list

   elseif nargout == 1,

      MODLIST = prune_modules(MODLIST);
      varargout{1} = MODLIST(end : -1 : 1);

   else

      error(msgid('Syntax:Error'), ...
           ['Call syntax is\n\t', ...
            mfilename, ' <command> [module list]  or\n\t', ...
            'mods = ', mfilename]);
   end
end

%--------------------------------------------------------------------------

function lst2 = prune_modules(lst)
   lst2 = {};

   if ~isempty(lst),
      pth  = split_path(path);
      mdir = module_dir;

      for i = 1 : numel(lst),
         d = fullfile(mdir, lst{i});

         assert (isdir(d) || isdir(lst{i}), ...
                 'Internal error defining module list.');

         if ~isdir(d), d = lst{i}; end

         if ~all(cellfun(@isempty, strfind(pth, d))),
            lst2 = [lst2, lst(i)];           %#ok  % Willfully ignore MLINT
         end
      end
   end
end

%--------------------------------------------------------------------------

function lst = add_modules(lst, mods)
   % Disallow mrstModule('add', ROOTDIR)
   %
   i = strcmp(ROOTDIR, mods) | cellfun(@isempty, mods);
   if any(i),
      dispif(mrstVerbose, ...
             ['mrstModule add ROOTDIR is explicitly ', ...
              'disallowed. Ignored.\n']);

      mods = mods(~i);
   end

   pth = [];
   for m = reshape(unique_modules(mods), 1, []),
      dirs = filter_module_dirs(m{1});

      if ~isempty(dirs),
         dirs = strcat(dirs, pathsep);
         pth  = [dirs{:}, pth];  %#ok

         i   = strcmpi(m{1}, lst);
         lst = [m, lst(~i)];
      end
   end

   if ~isempty(pth), addpath(pth); end
end

%--------------------------------------------------------------------------

function lst = clear_modules(lst)
   if ~isempty(lst),
      mroot = module_dir;
      for i = 1 : numel(lst),
         % Guard against false positives in module list.
         d = fullfile(mroot, lst{i});
         if isdir(d), lst{i} = d; end
      end

      pth  = split_path(path);
      patt = strcat('|', lst);
      patt = [ patt{:} ];
      patt = patt(2 : end);

      i = cellfun(@isempty, regexp(pth, patt));

      pth = strcat(pth(i), pathsep);
      path([ pth{:} ]);

      lst = {};
   end
end

%--------------------------------------------------------------------------

function dirs = filter_module_dirs(root)
   mroot = fullfile(module_dir, root);
   if isdir(mroot),
      root = mroot;
   elseif ~isdir(root),
      fprintf([root, ': No such module or directory!\n'])
      root = [];
   end

   if ~isempty(root),
      dirs  = split_path(genpath(root));
      vcdir = '\.(git|hg|svn)';

      is_vcdir = ~cellfun(@isempty, regexp(dirs, vcdir));
      exclude  = is_vcdir | cellfun(@isempty, dirs);

      dirs = dirs(~exclude);
   else
      dirs = {};
   end
end

%--------------------------------------------------------------------------

function d = module_dir()
   d = fullfile(ROOTDIR, 'modules');
end

%--------------------------------------------------------------------------

function pth = split_path(pth)
   try
      pth = regexp(pth, pathsep, 'split');
   catch  %#ok
      % Octave compatiblity.  It is an error to get here in an M run.
      pth = strsplit(pth, pathsep);
   end
end

%--------------------------------------------------------------------------

function mods = unique_modules(mods)
   excl_trail_fsep = @(s) s(1 : find(s ~= filesep, 1, 'last'));

   mods      = cellfun(excl_trail_fsep, ...
                       strtrim(mods)  , ...
                       'UniformOutput', false);
   [u, m, n] = unique(mods, 'last');  %#ok

   mods = mods(sort(m));
end
