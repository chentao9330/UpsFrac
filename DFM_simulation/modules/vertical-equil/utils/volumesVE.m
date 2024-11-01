function [totVol trappedVol freeVol] = volumesVE(G, sol, rock, fluid)
    % SYNOPSIS:
    %   [totVol trappedVol freeVol] = volumesVE(G, sol, rock, fluid)
    %   
    % PARAMETERS:
    % G         - 2D top surface grid used for VE-simulations
    % sol       - Solution state as defined by initResSolVE
    % rock      - rock for the top surface grid
    % fluid     - fluid object
    % 
    % RETURNS:
    % totVol, trappedVol, freeVol corresponding to total volume of CO2 in
    % the system, the trapped volume of CO2 and the free volume of CO2.
    freeVol = sum(sol.h.*rock.poro.*G.cells.volumes)*(1-fluid.sw);
    % Trapped CO2 is based on previous maximum of the CO2 column. If the
    % amount of CO2 has shrunk, some will have been trapped.
    trappedVol = sum((sol.h_max-sol.h).*rock.poro.*G.cells.volumes)*fluid.sr;
    totVol = trappedVol + freeVol;
end