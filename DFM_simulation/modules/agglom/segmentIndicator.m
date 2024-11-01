function [partition, edges] = segmentIndicator(G, indicator, numBins, varargin)
% Segments a fine grid into blocks according to indicator.
%
% SYNOPSIS:
%   partition = segmentIndicator(G, indicator, numBins)
%
% DESCRIPTION:
%   This function segments fine grid cells into the given number of bins
%   (numBins) according the given indicator value (indicator). The
%   groupings of cells are split into connected components.
%
% REQUIRED PARAMETERS:
%   G         - Grid data structure discretising the reservoir model
%               (fine grid, geological model).
%
%   indicator - Cell-wise value of some measure/indicator function.
%               Assuming positive numbers.
%
%   numBins   - Number of bins to segment the fine grid cells into.
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%
%   verbose  - Whether or not display number of blocks in the resulting
%              partition. Default value dependent upon global verbose
%              settings of function 'mrstVerbose'.
%
%
% RETURNS:
%   partition - Partition vector after segmenting all cells into blocks
%               according to the indicator
%
% SEE ALSO:
%   mergeBlocks, refineBlocks

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

opt  = struct('verbose' , mrstVerbose); % report number of blocks in segmentation?
opt = merge_options(opt, varargin{:});

edges = linspace(min(indicator), max(indicator), numBins + 1);
edges(end) = inf; % To avoid the maximum value of the indicator to end up in its own bin of number "numBins+1".

[N, partition] = histc(indicator, edges);

partition = compressPartition(partition);
partition = processPartition(G, partition, 'Reconnect', false);
partition = compressPartition(partition);

dispif(opt.verbose, 'SegmentIndicator: %d blocks\n', max(partition));

end
