% Supporting utilities for MATLAB Reservoir Simulation Toolbox (MRST).
%
% Files
%  ADI                   - ADI class: simple implementation of automatic differentiation for easy construction of jacobian matrices.
%  ROOTDIR               - Retrieve full path of Toolbox installation directory.
%  blockDiagIndex        - Compute subscript or linear index to nonzeros of block-diagonal matrix
%  buildmex              - Wrapper around MEX which abstracts away details of pathname generation.
%  cellFlux2faceFlux     - Transform cell-based flux field to face-based.
%  dinterpq1             - Compute derivative of piecewise linear interpolant.
%  dispif                - Produce textual output contingent upon predicate.
%  faceFlux2cellFlux     - Transform face-based flux field to cell-based.
%  faceFlux2cellVelocity - Transform face-based flux field to one constant velocity per cell.
%  findFilesSubfolders   - Search for files in a directory hierarchy
%  formatTimeRange       - Small utility which returns a human readable string from seconds.
%  geomspace             - Geometrically spaced vector.
%  initVariablesADI      - Initialize a set of automatic differentiation variables
%  invv                  - INVV(A, sz)
%  mcolon                - Compute vector of consecutive indices from separate lower/upper bounds.
%  md5sum                - md5sum - Compute md5 check sum of all input arguments
%  md5sum_fallback       - Alternative implementation of md5sum for systems without C compiler
%  merge_options         - Override default control options.
%  mrstDebug             - Globally control default settings for MRST debugging information.
%  mrstModule            - Query or modify list of activated add-on MRST modules
%  mrstVerbose           - Globally control default settings for MRST verbose information.
%  msgid                 - Construct Error/Warning message ID by prepending function name.
%  require               - Announce and enforce module dependency.
%  rldecode              - Decompress run length encoding of array A along dimension dim.
%  rlencode              - Compute run length encoding of array A along dimension dim.
%  ticif                 - Evaluate function TIC if input is true.
%  tocif                 - Evaluate function TOC if input is true.

%{
#COPYRIGHT#
%}
