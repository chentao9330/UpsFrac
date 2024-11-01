% Examples demonstrating the construction and manipulation of grid datastructures.
%
% Files
%  buildCornerPtNodes     - Construct physical nodal coordinates for CP grid.
%  buildCornerPtPillars   - Construct physical nodal coordinates for CP grid.
%  cart2active            - Compute active cell numbers from linear Cartesian index.
%  cartGrid               - Construct 2d or 3d Cartesian grid in physical space.
%  cellNodes              - Extract local-to-global vertex numbering for grid cells.
%  computeGeometry        - Compute geometry of grid.
%  extractSubgrid         - Construct valid grid definition from subset of existing grid cells.
%  grid_structure         - Grid structure used in MATLAB Reservoir Simulation Toolbox.
%  hexahedralGrid         - Construct valid grid definition from points and list of hexahedra
%  makeInternalBoundary   - Make internal boundary in grid along FACES
%  makeLayeredGrid        - Extrude 2D grid to layered 3D grid with n layers.
%  pebi                   - Compute dual grid of triangular grid G.
%  processFaults          - Construct fault structure from input specification (keyword 'FAULTS')
%  processGRDECL          - Compute grid topology and geometry from pillar grid description.
%  removeCells            - Remove cells from grid and renumber cells, faces and nodes.
%  removeFaultBdryFaces   - Remove fault faces on boundary
%  removeInternalBoundary - Remove internal boundary in grid by merging faces in face list N
%  removePinch            - Uniquify nodes, remove pinched faces and cells.
%  tensorGrid             - Construct Cartesian grid with variable physical cell sizes.
%  tetrahedralGrid        - Construct valid grid definition from points and tetrahedron list
%  triangleGrid           - Construct valid grid definition from points and triangle list

%{
#COPYRIGHT#
%}
