%% Relative permeability upscaling example
% This example demostrates upscaling of relative permeability on periodic
% grids. The example upscales a single block sampled from SPE10 using first
% permeability and then finds relative permeability based on viscous limit
% and various saturation values.
require spe10

%% Set up a simple grid with periodic boundaries
% Right -> Left
% Front -> Back
% Bottom -> Top

G   = cartGrid([5 5 2]);
G   = computeGeometry(G);

bcr{1}=pside([],G,'RIGHT',0);
bcr{2}=pside([],G,'FRONT',0);
bcr{3}=pside([],G,'BOTTOM',0);

bcl{1}=pside([],G,'LEFT',0);
bcl{2}=pside([],G,'BACK',0);
bcl{3}=pside([],G,'TOP',0);

dp = {0, 0, 0};

% Make periodic grid. We retain the regular grid for plotting, as plotGrid
% uses the boundary faces to plot grids: A fully periodic grid has, per
% definition, no boundary faces.

[Gp, bcp]=makePeriodicGridMulti3d(G, bcl, bcr, dp);
%% Extract a small subset of SPE10 to upscale.

x = 51; y = 11; z = 1;

rock = SPE10_rock(x:(x-1+G.cartDims(1)),...
                  y:(y-1+G.cartDims(2)),...
                  z:(z-1+G.cartDims(3)));
clf
plotCellData(G, log10(rock.perm(:,1)));
title('Fine scale permeability')

%% Do a single periodic upscale
% We upscale the permeability using two point flux approximation for the
% pressure solver
psolver = @(state, Grid, Fluid, BCP, Rock)...
           incompTPFA(state, Grid, computeTransGp(G, Grid, Rock),...
           Fluid, 'bcp', BCP);

% L is the size of the domain
L = max(G.faces.centroids) - min(G.faces.centroids);
% To find permeability we use a unitary fluid so that the mobility/relperm
% is equal to the saturation which is equal to one, removing any fluid
% specific effects.

fluid_pure = initSingleFluid('mu',1,'rho',1);

warning('off', 'mrst:periodic_bc')
perm2 = upscalePermeabilityPeriodic(Gp, bcp, 1, psolver, fluid_pure, rock, L);
warning('on', 'mrst:periodic_bc')

%% Load a two phase fluid for upscaling
% The data are synthetic and should not be used for anything but testing.
% rocklist contains a list of included property files in a simple format
% tabulated on water saturation.

current_dir = fileparts(mfilename('fullpath'));
fn    = fullfile(current_dir, 'rocklist.txt');

T  = readTabulatedJFluidFile(fn);

% Print the tabulated values from the first and only file
fprintf(' Sw         |Krw         |Kro         |J-func\n')
fprintf('--------------------------------------------------\n')
fprintf('%+1.4e |%+1.4e |%+1.4e |%+1.4e\n', T{1})

fluid = initSWOFFluidJfunc('mu' , [   10,  100] .* centi*poise     , ...
                             'rho', [1000, 600] .* kilogram/meter^3, ...
                             'table', T, ...
                             'satnum', 1, 'jfunc', true, 'rock', rock, ...
                             'surf_tens', 10*dyne/(centi*meter));
                         
%% Upscale relative permeability (viscous limit)
% We assume zero capillary forces and upscale using the viscous and capillary limit. 
[saturations_visc, kr_visc] = upscaleRelpermLimit(G, rock, fluid, 'type', 'fixed', 'limit', 'viscous');
[saturations_cap, kr_cap]   = upscaleRelpermLimit(G, rock, fluid, 'type', 'fixed', 'limit', 'capillary');

%% Plot the results for both limits
% The viscous limit is equal in all directions, while the capillary is not.
clf;
ph = {'water', 'oil'};
for i = 1:2
    subplot(2,1,i)
    hold on
    plot(saturations_visc, kr_visc{i});
    plot(saturations_cap, kr_cap{i}, '--.');
    title(['Relative permeability (Viscous/capillary limit), ' ph{i} ' phase']);
    xlabel('Saturation')
    legend({'x (viscous)', 'y (viscous)', 'z (viscous)'...
            'x (capillary)', 'y (capillary)', 'z (capillary)'}, 'location', 'West')
end

%% Set up a relative permeability upscaling run
% Saturations from 0:1 with resolution of ~20 data points.
saturations = 0:0.05:1;
% Pressure drop over periodic boundary used to induce flow.
dp_scale=1e-3;

% Ignore warnings from the implicit sovler as the solution is driven to
% steady state. It is natural that some steps fail during this process.
warning('off', 'implicitTransport:failure')
[sat_vec, kr, perm, krK] = upscaleRelperm(G, rock, fluid, dp_scale, saturations, 'periodic', false);
warning('on', 'implicitTransport:failure')

%% Plot the resulting relative permeability
% This is tabulated by water saturation in both cases
% As the default option is to use a pressure drop in x-direction, the
% x-values are significantly different from the y/z values which are
% similar, but not equal.

for i = 1:2
subplot(2,1,i)
plot(sat_vec, kr{i});
title(['Relative permeability, phase ' num2str(i)]);
xlabel('Water saturation')
legend({'x', 'y', 'z'}, 'location', 'West')
end

