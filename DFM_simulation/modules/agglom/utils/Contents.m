% Files
%  blockConnectivity            - Build block-to-neighbours map by transposing neighbourship definition
%  blockNeighbours              - Identify the neighbours of a particular coarse block.
%  blockNeighbourship           - Derive coarse-scale neighbourship from fine-scale information
%  coarse_sat                   - Converts a fine saturation field to a coarse saturation field, weighted
%  convertBC2Coarse             - Convert fine-scale boundary conditions to coarse scale.
%  convertRock2Coarse           - Create coarse-scale porosity field for solving transport equation.
%  convertSource2Coarse         - Accumulate fine-scale source terms to coarse scale
%  findConfinedBlocks           - Identify coarse blocks confined entirely within a single other block.
%  signOfFineFacesOnCoarseFaces - Identify fine-scale flux direction corresponding to coarse-scale outflow

%{
#COPYRIGHT#
%}
