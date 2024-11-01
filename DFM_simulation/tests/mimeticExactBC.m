%% z-direction
%  Check hydrostatic pressure solution in 5m water column with 
%  pressure = 0 to zero at z = 0
%
%  Check that both face pressures and cell pressures are exact.

G = cartGrid([9,10,11], [5,5,5]);
G = twister(G);
G = computeGeometry(G);
gravity off

rock.perm = ones(G.cells.num, 1)*darcy;
fluid     = initSingleFluid('mu' ,    1*centi*poise     , ...
                            'rho', 1000*kilogram/meter^3);
state     = initResSol(G, 0);
S     = computeMimeticIP(G, rock);
solver = @(bc) solveIncompFlow(state, G, S, fluid, 'bc', bc);


bc = [];
bc = pside(bc, G, 'bottom', 10);
bc = pside(bc, G, 'top',    0);
xr = solver(bc);
solution = @(z) 2*z;
mu       = fluid.properties(state);
flux     = darcy()*2/mu.*G.faces.normals(:,3);

n(1) = norm(xr.pressure     - solution(G.cells.centroids(:,3)), inf);
n(2) = norm(xr.facePressure - solution(G.faces.centroids(:,3)), inf);
n(3) = norm(xr.flux-flux);

fprintf('Max error in pressure, face pressure and face flux %d\n', max(n));



clf,hold on
title('Pressure as function of depth');
plot(G.faces.centroids(:,3), xr.facePressure, '*r');
plot(G.cells.centroids(:,3), xr.pressure, 'ob');
legend('face pressure','cell pressure','location','northwest');
z = linspace(0,5,100);
plot(z, solution(z), '-');
%axis([-0.1,5.1,-0.1,10.1])
hold off

if norm(n, inf) > 1, error('Wrong answer'); end





%% Plain copy for x-direction
bc = [];
bc = pside(bc, G, 'right', 10);
bc = pside(bc, G, 'left',    0);
xr = solver(bc);

solution = @(x) 2*x;
flux     = darcy()*2/mu.*G.faces.normals(:,1);

n(1) = norm(xr.pressure     - solution(G.cells.centroids(:,1)), inf);
n(2) = norm(xr.facePressure - solution(G.faces.centroids(:,1)), inf);
n(3) = norm(xr.flux-flux);

fprintf('Max error in pressure, face pressure and face flux %d\n', max(n));
if norm(n, inf) > 1, error('Wrong answer'); end



%% Plain copy for y-direction
bc = [];
bc = pside(bc, G, 'front', 10);
bc = pside(bc, G, 'back',    0);
xr = solver(bc);


solution = @(y) 2*y;
flux     = darcy()*2/mu.*G.faces.normals(:,1);

n(1) = norm(xr.pressure     - solution(G.cells.centroids(:,2)), inf);
n(2) = norm(xr.facePressure - solution(G.faces.centroids(:,2)), inf);
n(3) = norm(xr.flux-flux);

fprintf('Max error in pressure, face pressure and face flux %d\n', max(n));
if norm(n, inf) > 1, error('Wrong answer'); end
