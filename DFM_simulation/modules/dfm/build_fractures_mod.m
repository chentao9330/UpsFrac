function [vertices, edges, box] = build_fractures_mod (opts)
% Build fracture network that serves as constraints in the triangulation.
%
% Currently only odp-files are read, but other user-specified formats can
% easily be added here. 
%
% Copyright (C) 2006-2007 Uni Research AS
% This file is licensed under the GNU General Public License v3.0.


% filename is contained in the parameter; use the extension to
% determine the type of the file

filename = opts.fractures;
[~, ~, ext] = fileparts (filename);

% OpenOffice 2.0 presentations
ext = lower (ext(2:length (ext)));

[vertices, edges, box] = read_openoffice (filename, ext, opts);

% snap the fractures to the nearest grid; precision is relative to dimensions
vertices = snap_to_grid (vertices, box, opts);
