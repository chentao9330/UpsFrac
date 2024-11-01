function normS = twophaseTransport(doPlot, grav)
% twophaseTransport -- Examples of two-phase transport on small Cartesian
%                      grid with homogeneous, isotropic permeability.
%                      Compares explicit and implicit transport solvers.
%
% SYNOPSIS:
%           twophaseTransport()
%           twophaseTransport(doPlot)
%           twophaseTransport(doPlot, gravity)
%   normS = twophaseTransport(...)
%
% DESCRIPTION:
%   Function twophaseTransport runs the following examples:
%     1) Water flooding in a quarter five-spot well configuration on a
%        10-by-10 Cartesian grid with homogeneous, isotropic permeability
%        and porosity.  The reservoir is initially saturated with oil.
%
%     2) A gravity segregation example on the same 10-by-10 Cartesian grid
%        in which all cells initially have a water (and oil) saturation of
%        0.5.
%
%     3) A vertically separated gravity driven flow example on the same
%        10-by-10 Cartesian grid in which all cells to the left of the
%        separation line have an initial water saturation of 0.8, while all
%        cells to the right of the separation line have an initial water
%        saturation of 0.
%
% PARAMETERS:
%   doPlot  - Whether or not to plot results while the simulation is
%             running.  Logical.  Default value: doPlot = true.
%
%   gravity - Whether or not to consider gravity effects.
%             Examples 2) and 3) are only run when gravity is present.
%             Logical.  Default value: gravity = true.
%
% RETURNS:
%   normS - Euclidian norm of saturation differences, NORM(s_e - s_i), with
%           's_e' being the saturations computed with the explicit
%           transport solver and 's_i' being the saturations computed with
%           the implicit transport solver.  One scalar value for each test.

% $Date: 2011-02-08 23:20:03 +0100 (Tue, 08 Feb 2011) $
% $Revision: 6703 $

if nargin == 0,
   doPlot = true;
   grav   = true;

   G0 = gravity;
   gravity reset
end

verbose = mrstVerbose();


if verbose,
   if grav,
      disp('** Running tests of two-phase transport with gravity **');
   else
      disp('** Running test of two-phase transport **');
   end
   tic
end

%--------------------------------------------------------------------------
%- Initialize system ------------------------------------------------------
%
dims      = [10, 1, 10];
G         = cartGrid(dims);
G         = computeGeometry(G, 'Verbose', verbose);
rock.poro = repmat(0.3, [G.cells.num, 1]);
rock.perm = repmat(100*milli*darcy, [G.cells.num, 1]);
fluid     = initSimpleFluid('mu' , [   1,  10]*centi*poise     , ...
                            'rho', [1000, 700]*kilogram/meter^3, ...
                            'n'  , [   2,   2]);

gravity(grav);

S = computeMimeticIP(G, rock, 'Verbose', verbose);

%--------------------------------------------------------------------------
%- Run 3 cases if grav = true, 1 case otherwise ---------------------------
%
if grav, tests = 1:4; else tests = 1; end
normS = zeros([length(tests), 1]);

% Note on special MATLAB parsing quirk:
%   The MATLAB parser needs a little aid to unambiguously interpret the
%   body of an anonymous function if the definition is within a cell array.
%
%   Possible remedies are to:
%     - Remove all whitespace in the body of an anonymous function within a
%       cell array (which greatly impedes readability).
%     - Assign all anonymous functions to (separate) variables prior to
%       entering these variables into the cell array (which introduces a
%       great many "useless" variables).
%     - Enclose all anonymous functions in parentheses in order to separate
%       "readability whitespace" from essential, grouping-related
%       whitespace.  This is the choice made in the following definition of
%       variable 'fun'.
%
%   Further details may be found in the MATLAB documentation in the section
%      MATLAB -> Programming Fundamentals -> Types of Functions ->
%      Anonymous Functions -> Arrays of Anonymous Functions.
%
fun = {(@(xr, verb) test_bc  (xr, G, poreVolume(G,rock), verb)), ...
       (@(xr, verb) test_qfs (xr, G, poreVolume(G,rock), verb)), ...
       (@(xr, verb) test_grav(xr,                        verb)), ...
       (@(xr, verb) test_vss (xr, dims,                  verb))};

for i = tests,
   t = 0;
   [resSol, drive, tf, timestep] = fun{i}(initResSol(G, 0.0, [0,0]), verbose);

   %-----------------------------------------------------------------------
   %- Compute initial solution of pressure equation -----------------------
   %
   resSol      = solveIncompFlow(resSol, G, S, fluid, drive{:});
   resSol_impl = resSol;

   if doPlot,
      % Plot inital saturation
      figure;

      subplot(2,1,1),
         h = plotCellData(G, resSol.s(:,1));
         view([0, 0]);
         axis tight; axis equal;
         caxis([0, 1]);
         colorbar;
         
      subplot(2,1,2),
         h2 = plotCellData(G, resSol.s(:,1));
         view([0, 0]);
         axis tight; axis equal;
         caxis([0, 1]);
         colorbar;
         
   end

   %-----------------------------------------------------------------------
   %- Transport loop - solve saturation equation and update pressures -----
   %
   while t < tf,
      resSol      = explicitTransport(resSol, G, timestep, rock, fluid, ...
                                      drive{:});

      resSol_impl = implicitTransport(resSol_impl, G, timestep,     ...
                                      rock, fluid, drive{:},      ...
                                      'show_convergence',true);

      if max(max(abs([resSol.s(:,1), resSol_impl.s(:,1)]))) > 1 + 1.0e-5,
         disp('ERROR: ******* Saturation exceeds 1 **********')
         break
      end

      if doPlot,
         % Plot saturation
         delete(h);
         subplot(2,1,1)
            h = plotCellData(G, resSol.s(:,1));
            title(['Test ', int2str(i), ...
                   '. Water saturation - explicit solver'])
            view([0, 0]), axis tight equal
            caxis([0, 1]);

         delete(h2);
         subplot(2,1,2),
            h2 = plotCellData(G, resSol_impl.s(:,1));
            title('Water saturation - implicit solver')
            view([0, 0]), axis tight equal
            caxis([0 1]);
      end

      % Update solution of pressure equation.
      resSol      = solveIncompFlow(resSol     , G, S, fluid, drive{:});
      resSol_impl = solveIncompFlow(resSol_impl, G, S, fluid, drive{:});

      drawnow;
      t = t + timestep;
   end

   normS(i) = norm(resSol.s(:,1) - resSol_impl.s(:,1));
end

if grav && (exist('G0', 'var') == 1), gravity(G0); end

%--------------------------------------------------------------------------

function [resSol, drive, tf, dt] = test_bc(resSol, G, pv, verbose)
%- Quarter 5-spot --------------------------------------------
% Test with initial water saturation = 0 and injection of water.
%
if verbose,
   fprintf([ '\n\n****************************\n', ...
                 '** Test 1: bc             **\n', ...
                 '****************************\n']);
   tic
end
src   = [];
bc    = [];
%bc    = fluxside(bc, G, 'LEFT', 1e-5, 'sat', [1,0]);
bc    = pside(bc, G, 'LEFT',  51*barsa, 'sat', [1,0]);
bc    = pside(bc, G, 'RIGHT', 50*barsa, 'sat', [1,0]);
drive = {'src', src, 'bc', bc};

tf    = 10*convertFrom(sum(pv), day);
dt    = tf / 20;
%--------------------------------------------------------------------------

function [resSol, drive, tf, dt] = test_qfs(resSol, G, pv, verbose)
%- Quarter 5-spot --------------------------------------------
% Test with initial water saturation = 0 and injection of water.
%
if verbose,
   fprintf([ '\n\n****************************\n', ...
                 '** Test 2: Quarter 5-spot **\n', ...
                 '****************************\n']);
   tic
end
sourceCells = [ 1 ; G.cells.num ];
sources     = [ 1 ;    -1       ] .* meter^3 ./ day;  % m^3/s
src         = addSource([], sourceCells, sources, 'sat', [1,0; 1,0]);
bc          = [];
drive       = {'src', src, 'bc', bc};

tf          = convertFrom(sum(pv) / 2, day);
dt          = tf / 20;

%--------------------------------------------------------------------------

function [resSol, drive, tf, dt] = test_grav(resSol, verbose)
%- Gravity column ----------------------------------------
% Test of segregation flow caused by gravity.
% Initial water saturation in all cells:
%           ________
%          |        |
%          | s = 0.5|
%          |________|
%
if verbose,
   fprintf(['\n\n*********************************\n', ...
                '** Test 3: Gravity segregation **\n', ...
                '*********************************\n']);
   tic
end
src           = [];
bc            = [];
drive         = {'src', src, 'bc', bc};
resSol.s(:,1) = 0.5;

tf            = 3000 * day;
dt            = 100 * day;

%--------------------------------------------------------------------------

function [resSol, drive, tf, dt] = test_vss(resSol, dims, verbose)
%- Gravity vertical--------------------------------------
% Test with initial, vertically split saturation.
% Initial water saturation is:
%       ___________________
%      |          |        |
%      | s1 = 0.8 | s2 = 0 |
%      |__________|________|
%
if verbose,
   fprintf(['\n\n******************************************************\n', ...
                '** Test 4: Vertically separated gravity driven flow **\n', ...
                '******************************************************\n']);
   tic
end

src           = [];
bc            = [];
drive         = {'src', src, 'bc', bc};

resSol.s(:,1) = repmat(rldecode([0.8, 0], dims([1, 1]) ./ 2, 2) .', ...
                       [dims(3), 1]);
tf            = 3000 * day;
dt            = 100 * day;

