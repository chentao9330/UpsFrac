% This function calculates equivalent permeability based on files in a target directory.
% It is a refactored version of the CalcuKeq.m script.
%
% INPUT:
%   targetDir - The full path to the 'RESULT' folder containing the data files.
% Part of package: UpsFrac
% Author: Tao Chen et al.
% Email: chentao9330@gmail.com
% Copyright (c) 2024 UpsFrac Developers.
% All rights reserved.
% Updated: Oct 2024

%% --- User Settings ---
num_real = 1; % Total number of realizations to process
baseFolderName = 'reali_'; % Prefix for realization folder names
resultsSubFolder = 'RESULT'; % Subfolder name containing results
% ------------------

%% Define model parameters
Lx = 1000;          % Length of the model in the x-direction (m)
Ly = 1000;          % Length of the model in the y-direction (m)
dx = 100;           % Length of the coarse grid in the x-direction (m)
dy = 100;           % Length of the coarse grid in the y-direction (m)
nx = Lx / dx;
ny = Ly / dy;
perm_matrix_x = 9.87e-16;
perm_matrix_y = 9.87e-16;

% Get current main script directory
baseDir = pwd;

% Start timer
tic;

% Create results summary folder if it doesn't exist
resuSumFolder = fullfile(baseDir, 'resu_sum');
if ~isfolder(resuSumFolder)
    mkdir(resuSumFolder);
    fprintf('Created results summary folder: %s\n', resuSumFolder);
end

fprintf('Starting batch processing...\n');

% Loop through each realization
for i = 1:num_real
    % Build full path to target 'RESULT' folder
    targetDir = fullfile(baseDir, [baseFolderName, num2str(i)], resultsSubFolder);
    
    % Check if target folder exists
    if isfolder(targetDir)
        fprintf('\nProcessing realization %d in folder: %s\n', i, targetDir);
        
        % --- Core function call ---
        try
            CalcuKeq_function(targetDir, Lx, Ly, dx, dy, nx, ny, perm_matrix_x, perm_matrix_y);
            
            % --- Copy and rename epcomf.txt to resu_sum folder ---
            sourceFile = fullfile(targetDir, 'epcomf.txt');
            destFile = fullfile(resuSumFolder, sprintf('epcomf_%s%d.txt', baseFolderName, i));
            
            if exist(sourceFile, 'file') == 2
                copyfile(sourceFile, destFile);
                fprintf('Copied results to: %s\n', destFile);
            else
                fprintf('Warning: epcomf.txt not found in %s\n', targetDir);
            end
            
        catch ME
            fprintf('Error occurred in realization %d: %s\n', i, ME.message);
            % Continue processing other realizations
            continue;
        end
        
    else
        fprintf('\nWarning: Directory not found, skipping: %s\n', targetDir);
    end
end

fprintf('\nAll realizations processed successfully!\n');
toc; % Display total elapsed time




function CalcuKeq_function(targetDir, Lx, Ly, dx, dy, nx, ny, perm_matrix_x, perm_matrix_y)
% This function calculates equivalent permeability based on files in a target directory.
% It is a refactored version of the CalcuKeq.m script.
%
% INPUT:
%   targetDir - The full path to the 'RESULT' folder containing the data files.
%   Lx, Ly - Length of the model in x and y directions (m)
%   dx, dy - Length of the coarse grid in x and y directions (m)
%   nx, ny - Number of grids in x and y directions
%   perm_matrix_x, perm_matrix_y - Matrix permeability values

% Use the model parameters passed as function arguments

format long e;
% Note: Variable initializations from the original script are kept,
% although many are re-assigned within the loop.

%% Open main output files using full paths
voidgrid_path = fullfile(targetDir, 'voidgrid.txt');
epcomf_path = fullfile(targetDir, 'epcomf.txt');
fid10000 = fopen(voidgrid_path, 'w');
fid20000 = fopen(epcomf_path, 'w');

% This variable holds the full path to the output file for convenience
comt = epcomf_path;

%% Loop through all grid numbers
for num_grid = 1:nx * ny
    
    seq = sprintf('%d', num_grid);
    post = '.txt';
    pre_f = 'f.txt';
    
    % Construct full paths for all potential input files
    cf = fullfile(targetDir, strcat(strtrim(seq), pre_f));
    hy = fullfile(targetDir, strcat('hybrid', strtrim(seq), post));
    cen = fullfile(targetDir, strcat('centroids', strtrim(seq), post));
    flua = fullfile(targetDir, strcat('ffluxpro1', strtrim(seq), post));
    flub = fullfile(targetDir, strcat('ffluxpro2', strtrim(seq), post));
    
    % Set default matrix permeability for this grid
    kxx = perm_matrix_x;
    kxy = 0;
    kyx = 0;
    kyy = perm_matrix_y;
    
    % Check if the primary fracture data file exists for this grid
    if exist(flua, 'file') == 2
        
        % Assign full-path variables for clarity within the calculation block
        frac = cf;
        hybrid = hy;
        centroids = cen;
        ffluxpro1 = flua;
        ffluxpro2 = flub;
        comftxt = comt; % This is already a full path
        matrixp = kxx;
        
        % Initialize variables for reading files
        fracn = 0;
        fluxn = 0;
        row = 7;
        gn = 4;
        
        % --- Pre-read to get dimensions ---
        % Read number of fractures (lines in frac file)
        fid = fopen(frac, 'r');
        while ~feof(fid)
            buffer = fgetl(fid);
            if ischar(buffer)
                fracn = fracn + 1;
            end
        end
        fclose(fid);
        
        % Read number of fluxes (lines in centroids file)
        fid = fopen(centroids, 'r');
        while ~feof(fid)
            puffer = fgetl(fid);
            if ischar(puffer)
                fluxn = fluxn + 1;
            end
        end
        fclose(fid);
        
        % Pre-allocate memory for speed
        g = zeros(gn, fracn);
        f = zeros(row, fracn);
        
        % --- Process fracture geometry data ---
        % Create a temporary file 'fra.txt' in the target directory
        fra_txt_path = fullfile(targetDir, 'fra.txt');
        fid = fopen(frac, 'r');
        fid2 = fopen(fra_txt_path, 'w');
        
        data = fscanf(fid, '%f %f %f %f %f', [5 inf]);
        fprintf(fid2, '%f %f %f %f\n', data(1:4, :)); % Write only first 4 columns
        
        fclose(fid);
        fclose(fid2);
        
        % Read the processed fracture data
        fid = fopen(fra_txt_path, 'r');
        for i = 1:fracn
            data = fscanf(fid, '%f %f %f %f', [4 1]);
            f(1:4, i) = data;
            f(5, i) = (f(4, i) - f(2, i)) / (f(3, i) - f(1, i)); % Slope
            f(6, i) = cos(atan(f(5, i))); % cos(angle)
            f(7, i) = sin(atan(f(5, i))); % sin(angle)
        end
        fclose(fid);
        
        % --- Main calculation loop for two pressure drop problems (x and y) ---
        for p = 1:2
            prob = p;
            
            if prob == 1
                % Problem 1: Pressure drop in X direction
                
                % Allocate memory for structures
                upsca = repmat(struct('fn', 0, 'bound', '', 'flux', zeros(1, 7)), fluxn, 1);
                
                % Initialize flux accumulators
                tqmx = 0; tqfx = 0;
                tqmy = 0; tqfy = 0;
                
                % Read hybrid, centroids, and ffluxpro1 files
                fid1 = fopen(hybrid, 'r');
                fid2 = fopen(centroids, 'r');
                fid3 = fopen(ffluxpro1, 'r');
                for i = 1:fluxn
                    upsca(i).fn = i;
                    upsca(i).flux(1) = fscanf(fid1, '%f', 1);       % is_frac
                    upsca(i).flux(2:3) = fscanf(fid2, '%f %f', 2);  % centroid x, y
                    upsca(i).flux(4) = fscanf(fid3, '%e', 1);       % flux value
                end
                fclose(fid1); fclose(fid2); fclose(fid3);
                
                % Match fracture geometry to boundary intersection points
                for i = 1:fracn
                    ia = f(1, i); ib = f(2, i); ic = f(3, i); id = f(4, i);
                    for j = 1:fluxn
                        if upsca(j).flux(1) == 1 % If it's a fracture boundary
                            fa = upsca(j).flux(2); fb = upsca(j).flux(3);
                            if (abs(fa - ia) < 0.005 && abs(fb - ib) < 0.005) || ...
                               (abs(fa - ic) < 0.005 && abs(fb - id) < 0.005)
                                upsca(j).flux(5:7) = f(5:7, i); % Copy slope and orientation
                            end
                        end
                    end
                end
                
                % Calculate total flux across boundaries
                for i = 1:fluxn
                    if upsca(i).flux(4) > 0 % Only consider outgoing flux
                        if upsca(i).flux(1) == 0 % Matrix flux
                            if upsca(i).flux(2) == dx
                                tqmx = tqmx + upsca(i).flux(4);
                            elseif upsca(i).flux(3) == 0
                                tqmy = tqmy - upsca(i).flux(4);
                            elseif upsca(i).flux(3) == dy
                                tqmy = tqmy + upsca(i).flux(4);
                            end
                        else % Fracture flux
                            if upsca(i).flux(2) == dx
                                tqfx = tqfx + upsca(i).flux(4) * upsca(i).flux(6); % q*cos
                                tqfy = tqfy + upsca(i).flux(4) * upsca(i).flux(7); % q*sin
                            elseif upsca(i).flux(3) == 0
                                if upsca(i).flux(5) > 0
                                    tqfx = tqfx - upsca(i).flux(4) * upsca(i).flux(6);
                                else
                                    tqfx = tqfx + upsca(i).flux(4) * upsca(i).flux(6);
                                end
                                tqfy = tqfy - abs(upsca(i).flux(7)) * upsca(i).flux(4);
                            elseif upsca(i).flux(3) == dy
                                if upsca(i).flux(5) > 0
                                    tqfx = tqfx + upsca(i).flux(4) * upsca(i).flux(6);
                                else
                                    tqfx = tqfx - upsca(i).flux(4) * upsca(i).flux(6);
                                end
                                tqfy = tqfy + abs(upsca(i).flux(7)) * upsca(i).flux(4);
                            end
                        end
                    end
                end
                
                qx = tqfx + tqmx;
                qy = tqfy + tqmy;
                
                kxx = (qx * 0.001) / dx;
                kyx = (qy * 0.001) / dx;
                
                clear upsca;
                
            elseif prob == 2
                % Problem 2: Pressure drop in Y direction
                
                upscb = repmat(struct('fn', 0, 'bound', '', 'flux', zeros(1, 7)), fluxn, 1);
                
                % Initialize flux accumulators
                tqmx = 0; tqfx = 0;
                tqmy = 0; tqfy = 0;
                
                % Read hybrid, centroids, and ffluxpro2 files
                fid1 = fopen(hybrid, 'r');
                fid2 = fopen(centroids, 'r');
                fid4 = fopen(ffluxpro2, 'r');
                for i = 1:fluxn
                    upscb(i).fn = i;
                    upscb(i).flux(1) = fscanf(fid1, '%f', 1);
                    upscb(i).flux(2:3) = fscanf(fid2, '%f %f', 2);
                    upscb(i).flux(4) = fscanf(fid4, '%e', 1);
                end
                fclose(fid1); fclose(fid2); fclose(fid4);
                
                % Match fracture geometry to boundary intersection points
                for i = 1:fracn
                    ia = f(1, i); ib = f(2, i); ic = f(3, i); id = f(4, i);
                    for j = 1:fluxn
                        if upscb(j).flux(1) == 1
                            fa = upscb(j).flux(2); fb = upscb(j).flux(3);
                            if (abs(fa - ia) < 0.005 && abs(fb - ib) < 0.005) || ...
                               (abs(fa - ic) < 0.005 && abs(fb - id) < 0.005)
                                upscb(j).flux(5:7) = f(5:7, i);
                            end
                        end
                    end
                end
                
                % Calculate total flux across boundaries
                for i = 1:fluxn
                    if upscb(i).flux(4) > 0 % Only consider outgoing flux
                        if upscb(i).flux(1) == 0 % Matrix flux
                            if upscb(i).flux(2) == 0
                                tqmx = tqmx - upscb(i).flux(4);
                            elseif upscb(i).flux(2) == dx
                                tqmx = tqmx + upscb(i).flux(4);
                            elseif upscb(i).flux(3) == dy
                                tqmy = tqmy + upscb(i).flux(4);
                            end
                        else % Fracture flux
                            if upscb(i).flux(2) == 0
                                tqfx = tqfx - upscb(i).flux(4) * upscb(i).flux(6);
                                tqfy = tqfy - upscb(i).flux(4) * upscb(i).flux(7);
                            elseif upscb(i).flux(2) == dx
                                tqfx = tqfx + upscb(i).flux(4) * upscb(i).flux(6);
                                tqfy = tqfy + upscb(i).flux(4) * upscb(i).flux(7);
                            elseif upscb(i).flux(3) == dy
                                if upscb(i).flux(5) > 0
                                    tqfx = tqfx + upscb(i).flux(4) * upscb(i).flux(6);
                                else
                                    tqfx = tqfx - upscb(i).flux(4) * upscb(i).flux(6);
                                end
                                tqfy = tqfy + abs(upscb(i).flux(7)) * upscb(i).flux(4);
                            end
                        end
                    end
                end
                
                qx = tqfx + tqmx;
                qy = tqfy + tqmy;
                
                kxy = (qx * 0.001) / dx;
                kyy = (qy * 0.001) / dx;
                
                clear upscb;
                
                % --- Final Sanity Checks and Write Results ---
                if kxx <= 0 || kxx > 1
                    kxx = matrixp;
                end
                if kyy <= 0 || kyy > 1
                    kyy = matrixp;
                end
                if abs(kxy) < matrixp || abs(kxy) > 1
                    kxy = 0;
                end
                if abs(kyx) < matrixp || abs(kyx) > 1
                    kyx = 0;
                end
                
                [~, fname, fext] = fileparts(frac);
                fprintf(fid20000, '%s %0.8e %0.8e %0.8e %0.8e\n', [fname, fext], kxx, kxy, kyx, kyy);
            end
        end
        
    else
        % This block executes if the grid has no fractures (flua file does not exist)
        [~, fname, fext] = fileparts(cf);
        fprintf(fid20000, '%s %0.8e %0.8e %0.8e %0.8e\n', [fname, fext], kxx, kxy, kyx, kyy);
        fprintf(fid10000, 'grid %d does not exist.\n', num_grid);
    end
end

%% Close files and clean up
fclose(fid10000);
fclose(fid20000);

% Delete the temporary file
if exist(fra_txt_path, 'file')
    delete(fra_txt_path);
end

% Display a completion message in the command window
fprintf('Processing finished for directory: %s\n', targetDir);

end