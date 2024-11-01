%%
%  Check hydrostatic pressure solution in 5m water column with 
%  pressure = 0 to zero at z = 0
%
%  Check that both face pressures and cell pressures are exact.

G = cartGrid([9,10,11], [5,5,5]);
G = twister(G);
G = computeGeometry(G);

grav = gravity;
gravity([0, 0, 10]);
gravity on

bc = [];
bc = pside(bc, G, 'top', 0);

rock.perm = ones(G.cells.num, 1)*milli*darcy();
fluid     = initSingleFluid('mu' ,    1*centi*poise     , ...
                            'rho', 1000*kilogram/meter^3);
state     = initResSol(G, 0);
T         = computeTrans(G, rock);

xr    = incompTPFA(state, G, T, fluid, 'bc', bc);
solution = @(z) 1e4*z;

n(1) = norm(xr.pressure     - solution(G.cells.centroids(:,3)));
n(2) = norm(xr.facePressure - solution(G.faces.centroids(:,3)));
n(3) = norm(xr.flux);

fprintf('Max error in pressure, face pressure and face flux %d\n', max(n));


clf,hold on
title('Pressure as function of depth');
plot(G.faces.centroids(:,3), xr.facePressure, '*r');
plot(G.cells.centroids(:,3), xr.pressure, 'ob');
legend('face pressure','cell pressure','location','northwest');
z = linspace(0,5,100);
plot(z, solution(z), '-');
axis([-0.1,5.1,-1e4,5.1e4])
hold off

if norm(n, inf) > 1, error('Wrong answer'); end

gravity(grav);
