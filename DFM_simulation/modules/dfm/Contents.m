% Discrete Fracture Matrix (DFM) module
%
% Routines supporting the DFM method. Examples in /tests
% Add module by typing MrstModule add "path of DFM folder"
% Tor Harald Sandve 
%
% new files
%   addhybrid              - Generates the hybrid cells
%   computeHybridTrans     - Computes the hybrid transmissibilities using TPFA 
%   computeMultiPointTrans - Computes the hybrid transmissibilities using
%                            O-MPFA
%   testNormals            - Rutine to check that the direction of the normals
%                            is from neighbor 1 to 2. 
%   /plotting/plotFractures -Plots 2d hybridcells. 
%
% Files modified from  core MRST functions. Use these to 
%
%   computeTrans_DFM           - Modified to compute hybrid-normal
%                                transmissibilities using a modified TPFA
%   computeMultiPointTrans_DFM - Modified to compute hybrid-normal
%                                transmissibilities using a modified MPFA 
%   computeTPFA_DFM            - Modified to allow for cell2cell connections
%   computeMPFA_DFM            - Modified to allow for cell2cell connections
%   twophaseJacobian_DFM       - Modified to allow for cell2cell connections
%   explicitTransport_DFM      - Modified to allow for cell2cell connections
%   implicitTransport_DFM      - Modified to allow for cell2cell connections
%   removeInternalBoundary_DFM - Modified to update grid geometry as well as the topology  
%   /plotting/                 - Modification to plotting rutines to handle
%                                the hybrid cells. 
%
%
% /localUnmodifiedFiles. 
%  - Files needed since we can not acces private folders from this location
%   
% /examples:
%   simpleFracturedMedia2D - A simple fracture system in a Cartesian 2D
%                            grid. The extension to 3D is straightforward
% 
