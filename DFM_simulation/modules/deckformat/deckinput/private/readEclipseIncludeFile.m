function varargout = readEclipseIncludeFile(fun, fid, dirname, varargin)
%Read an ECLIPSE INCLUDE file.
%
% SYNOPSIS:
%   [ret{1:nret}] = readEclipseIncludeFile(fun, fid, dirname, ...)
%
% PARAMETERS:
%   fun     - Callback function handle.  Assumed to support the syntax
%
%                [ret{1:nret}] = fun(fid, dirname, ...)
%
%             Function 'fun' will be called with 'fid' set to the FOPEN
%             return value of the INCLUDE file name, while 'dirname' will
%             be the complete directory name of the INCLUDE file name.
%
%   fid     - Valid file identifier as obtained by FOPEN.  The file pointer
%             FTELL(fid) is assumed to be placed directly after the
%             'INCLUDE' keyword and strictly before the name of the file
%             which will be INCLUDEd.
%
%             Function 'readEclipseIncludeFile' will, upon successful
%             return, place FTELL(fid) directly after the keyword-closing
%             slash character.
%
%   dirname - Complete directory name of file from which the input file
%             identifier 'fid' was derived through FOPEN.
%
%   ...     - Additional function parameters.  These parameters will be
%             passed unchanged on to function 'fun'.
%
% RETURNS:
%   ret     - Any and all return values from function 'fun'.  It is the
%             responsibility of the caller of 'readEclipseIncludeFile' to
%             supply sufficient number of output arrays (i.e., 'nret') to
%             store these return values.

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

   lin = '';
   while ischar(lin) && isempty(strtrim(lin)),
      lin = fgetl(fid);
      lin = regexprep(lin, '--.*$', '');
   end

   if ischar(lin),
      lin = strtrim(lin);

      if any(lin == ''''),
         if sum(lin == '''') ~= 2,
            fclose(fid);
            error('INCLUDE argument not correctly delimited.');
         end

         % Extract path-name portion of INCLUDE argument.
         inc_fn = regexprep(lin, '.*''([^'']+)''.*', '$1');
      else
         inc_fn = sscanf(lin, '%s');
      end

      p = find(lin == '/', 1, 'last');
      if isempty(p),
         terminated = false;
      else
         terminated = p == numel(lin) || isspace(lin(p + 1));
      end
   end

   if inc_fn(end) == '/', inc_fn = inc_fn(1 : end - 1); end

   % Gobble up keyword-closing '/' character if not already read.
   if ~terminated,
      p     = ftell(fid);
      fn    = fopen(fid);

      slash = fscanf(fid, '%s', 1);  % Possibly too weak.

      if ~strcmp(slash, '/'),
         fclose(fid);
         error(msgid('Include:WrongfulTermination'), ...
              ['INCLUDE keyword not correctly terminated at ', ...
               'position %lu in file ''%s'''], p, fn);
      end
   end

   inc_fn(inc_fn == '/') = filesep;  % FILESEP is always a single char.
   if inc_fn(1) ~= filesep,
      % Translate relative pathname to absolute pathname.
      inc_fn = fullfile(dirname, inc_fn);
   end

   [inc_fid, msg] = fopen(inc_fn, 'rt');
   if inc_fid < 0, error([inc_fn, ': ', msg]); end

   try
      % Call back to our (likely) caller with the new file.
      [varargout{1:nargout}] = fun(inc_fid, ...
                                   dirname, ...  % or fileparts(inc_fn)
                                   varargin{:});
   catch %#ok
      err = lasterror;  %#ok
      try  %#ok
         % Don't leak fids, but don't gripe about a child already closing
         % the fid.
         fclose(inc_fid);
      end
      rethrow(err);
   end
   fclose(inc_fid);
end
