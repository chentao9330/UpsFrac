% Part of package: UpsFrac
% Author: Tao Chen et al.
% Email: chentao9330@gmail.com
% Copyright (c) 2024 UpsFrac Developers.
% All rights reserved.
% Updated: Oct 2024




clear all;
clc;

%%
% Set the values for lx, ly, dx, dy, and matrix permeability
% Define the parameters for the model
lx = 1000;        % Length of the model in the x-direction (m)
ly = 1000;        % Length of the model in the y-direction (m)
dx = 100;         % Length of the coarse grid in the x-direction (m)
dy = 100;         % Length of the coarse grid in the y-direction (m)
nx=lx/dx;
ny=ly/dy;
perm_matrix_x=9.87e-16;
perm_matrix_y=9.87e-16;



format long e;
% Initialize variables
cf = '';
hy = '';
cen = '';
flua = '';
flub = '';
pre_f = '';
pre_h = '';
pre_c = '';
pre_fla = '';
pre_flb = '';
seq = '';
post = '';
comt = '';
num_grid = 0;
alive = false;
kxx = 0;
kxy = 0;
kyx = 0;
kyy = 0;

% Open files
fid10000 = fopen('voidgrid.txt', 'w');
fid20000 = fopen('epcomf.txt', 'w');

comt = 'epcomf.txt';

% Loop through grid numbers
for num_grid = 1:nx*ny
    
    
    seq = sprintf('%d', num_grid);
    post = '.txt';
    pre_f = 'f.txt';
    
    cf = strcat(strtrim(seq), pre_f);
    hy = strcat('hybrid', strtrim(seq), post);
    cen = strcat('centroids', strtrim(seq), post);
    flua =strcat( 'ffluxpro1' , strtrim(seq), post);
    flub = strcat('ffluxpro2', strtrim(seq), post);
    
    kxx = perm_matrix_x;
    kxy = 0;
    kyx = 0;
    kyy = perm_matrix_y;
    
    % Check if file exists
    if exist(flua, 'file') == 2
        
        frac=cf;hybrid=hy;centroids=cen;ffluxpro1=flua;ffluxpro2=flub;comftxt=comt;matrixp=kxx;
        
        % Initialize variables
        fracn = 0;
        fluxn = 0;
        row = 7;
        gn = 4;
        dx = dx;
        dy = dy;
        g = zeros(gn, fracn);
        f = zeros(row, fracn);
        upsca = struct('fn', {}, 'bound', {}, 'flux', {});
        upscb = struct('fn', {}, 'bound', {}, 'flux', {});
        
        % Read frac file
        fid = fopen(frac, 'r');
        while ~feof(fid)
            buffer = fgetl(fid);
            if ischar(buffer)
                fracn = fracn + 1;
            end
        end
        fclose(fid);
        
        % Read centroids file
        fid = fopen(centroids, 'r');
        while ~feof(fid)
            puffer = fgetl(fid);
            if ischar(puffer)
                fluxn = fluxn + 1;
            end
        end
        fclose(fid);
        
        % Allocate memory
        g = zeros(gn, fracn);
        f = zeros(row, fracn);
        
        
        % Read frac file again
        fid = fopen(frac, 'r');
        fid2 = fopen('fra.txt', 'w');
        
        data = fscanf(fid, '%f %f %f %f %f', [5 inf]);
        fprintf(fid2, '%f %f %f %f\n', data);
        
        fclose(fid);
        fclose(fid2);
        
        
        % Read fra.txt file
        fid = fopen('fra.txt', 'r');
        for i = 1:fracn
            data = fscanf(fid, '%f %f %f %f', [4 1]);
            %disp(data);
            f(1:4, i) = data;
            f(5, i) = (f(4, i) - f(2, i)) / (f(3, i) - f(1, i));
            f(6, i) = cos(atan(f(5, i)));
            f(7, i) = sin(atan(f(5, i)));
        end
        fclose(fid);
        
        for p = 1:2
            prob = p;
            
            
            if prob == 1
                
                % Allocate memory for upsca and upscb
                upsca = repmat(struct('fn', 0, 'bound', '', 'flux', zeros(1, 7)), fluxn, 1);
                upscb = repmat(struct('fn', 0, 'bound', '', 'flux', zeros(1, 7)), fluxn, 1);
                
                % Initialize variables
                qmx = 0;
                tqmx = 0;
                qfx = 0;
                tqfx = 0;
                qx = 0;
                qmy = 0;
                tqmy = 0;
                qfy = 0;
                tqfy = 0;
                qy = 0;
                
                % Read hybrid, centroids, and ffluxpro1 files
                fid1 = fopen(hybrid, 'r');
                fid2 = fopen(centroids, 'r');
                fid3 = fopen(ffluxpro1, 'r');
                for i = 1:fluxn
                    upsca(i).fn = i;
                    upsca(i).flux(1) = fscanf(fid1, '%f', 1);
                    upsca(i).flux(2:3) = fscanf(fid2, '%f %f', 2);
                    upsca(i).flux(4) = fscanf(fid3, '%e', 1);
                end
                fclose(fid1);
                fclose(fid2);
                fclose(fid3);
                
                
                % Add orientation information to fracture boundaries
                for i = 1:fracn
                    ia = f(1, i);
                    ib = f(2, i);
                    ic = f(3, i);
                    id = f(4, i);
                    for j = 1:fluxn
                        if upsca(j).flux(1) == 1
                            
                            fa = upsca(j).flux(2);
                            fb = upsca(j).flux(3);
                            if abs(fa - ia) < 0.005 && abs(fb - ib) < 0.005
                                upsca(j).flux(5:7) = f(5:7, i);
                            elseif abs(fa - ic) < 0.005 && abs(fb - id) < 0.005
                                upsca(j).flux(5:7) = f(5:7, i);
                            end
                            
                        end
                    end
                end
                
                % Output the flux on boundaries
                
                for i = 1:fluxn
                    if upsca(i).flux(4) > 0
                        if upsca(i).flux(1) == 0
                            if upsca(i).flux(2) == 0
                                upsca(i).bound = 'left';
                            elseif upsca(i).flux(2) == dx
                                upsca(i).bound = 'right';
                                qmx = upsca(i).flux(4);
                                tqmx = tqmx + qmx;
                            elseif upsca(i).flux(3) == 0
                                upsca(i).bound = 'back';
                                qmy = -upsca(i).flux(4);
                                tqmy = tqmy + qmy;
                            elseif upsca(i).flux(3) == dy
                                upsca(i).bound = 'front';
                                qmy = upsca(i).flux(4);
                                tqmy = tqmy + qmy;
                            end
                        else % fracture flux
                            if upsca(i).flux(2) == 0
                                upsca(i).bound = 'left';
                                
                            elseif upsca(i).flux(2) == dx
                                upsca(i).bound = 'right';
                                qfx = upsca(i).flux(4) * upsca(i).flux(6);
                                tqfx = tqfx + qfx;
                                qfy = upsca(i).flux(4) * upsca(i).flux(7);
                                tqfy = tqfy + qfy;
                                %                                 disp('right');
                                
                            elseif upsca(i).flux(3) == 0
                                upsca(i).bound = 'back';
                                if upsca(i).flux(5) > 0
                                    qfx = -upsca(i).flux(4) * upsca(i).flux(6);
                                    tqfx = tqfx + qfx;
                                elseif upsca(i).flux(5) < 0
                                    qfx = upsca(i).flux(4) * upsca(i).flux(6);
                                    tqfx = tqfx + qfx;
                                end
                                qfy = -abs(upsca(i).flux(7)) * upsca(i).flux(4);
                                tqfy = tqfy + qfy;
                                %                                 disp('back');
                            elseif upsca(i).flux(3) == dy
                                upsca(i).bound = 'front';
                                if upsca(i).flux(5) > 0
                                    qfx = upsca(i).flux(4) * upsca(i).flux(6);
                                    tqfx = tqfx + qfx;
                                elseif upsca(i).flux(5) < 0
                                    qfx = -upsca(i).flux(4) * upsca(i).flux(6);
                                    tqfx = tqfx + qfx;
                                end
                                %                                 disp('front');
                                %                                 disp(upsca(i));
                                qfy = abs(upsca(i).flux(7)) * upsca(i).flux(4);
                                tqfy = tqfy + qfy;
                                
                                
                            end
                        end
                    end
                end
                
                
                qx = tqfx + tqmx;
                qy = tqfy + tqmy;
                
                kxx = (qx * 0.001) / dx;
                kyx = (qy * 0.001) / dx;
                
                %                Deallocate memory
                clear upsca upscb;
                
                qx = 0;
                qy = 0;
                
                
            elseif prob == 2
                
                
                % Allocate memory for upsca and upscb
                upsca = repmat(struct('fn', 0, 'bound', '', 'flux', zeros(1, 7)), fluxn, 1);
                upscb = repmat(struct('fn', 0, 'bound', '', 'flux', zeros(1, 7)), fluxn, 1);
                
                % Initialize variables
                qmx = 0;
                tqmx = 0;
                qfx = 0;
                tqfx = 0;
                qx = 0;
                qmy = 0;
                tqmy = 0;
                qfy = 0;
                tqfy = 0;
                qy = 0;
                
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
                fclose(fid1);
                fclose(fid2);
                fclose(fid4);
                
                % Add orientation information to fracture boundaries
                for i = 1:fracn
                    ia = f(1, i);
                    ib = f(2, i);
                    ic = f(3, i);
                    id = f(4, i);
                    for j = 1:fluxn
                        if upscb(j).flux(1) == 1
                            fa = upscb(j).flux(2);
                            fb = upscb(j).flux(3);
                            if abs(fa - ia) < 0.005 && abs(fb - ib) < 0.005
                                upscb(j).flux(5:7) = f(5:7, i);
                            elseif abs(fa - ic) < 0.005 && abs(fb - id) < 0.005
                                upscb(j).flux(5:7) = f(5:7, i);
                            end
                        end
                    end
                end
                
                % Output the flux on boundaries
                fid5 = fopen('fluxm.txt', 'w');
                fid7 = fopen('fluxf.txt', 'w');
                for i = 1:fluxn
                    if upscb(i).flux(4) > 0
                        if upscb(i).flux(1) == 0
                            if upscb(i).flux(2) == 0
                                upscb(i).bound = 'left';
                                qmx = -upscb(i).flux(4);
                                tqmx = tqmx + qmx;
                            elseif upscb(i).flux(2) == dx
                                upscb(i).bound = 'right';
                                qmx = upscb(i).flux(4);
                                tqmx = tqmx + qmx;
                            elseif upscb(i).flux(3) == 0
                                upscb(i).bound = 'back';
                            elseif upscb(i).flux(3) == dy
                                upscb(i).bound = 'front';
                                qmy = upscb(i).flux(4);
                                tqmy = tqmy + qmy;
                            end
                        else
                            if upscb(i).flux(2) == 0
                                upscb(i).bound = 'left';
                                qfx = -upscb(i).flux(4) * upscb(i).flux(6);
                                tqfx = tqfx + qfx;
                                qfy = -upscb(i).flux(4) * upscb(i).flux(7);
                                tqfy = tqfy + qfy;
                            elseif upscb(i).flux(2) == dx
                                upscb(i).bound = 'right';
                                qfx = upscb(i).flux(4) * upscb(i).flux(6);
                                tqfx = tqfx + qfx;
                                qfy = upscb(i).flux(4) * upscb(i).flux(7);
                                tqfy = tqfy + qfy;
                            elseif upscb(i).flux(3) == 0
                                upscb(i).bound = 'back';
                            elseif upscb(i).flux(3) == dy
                                upscb(i).bound = 'front';
                                if upscb(i).flux(5) > 0
                                    qfx = upscb(i).flux(4) * upscb(i).flux(6);
                                    tqfx = tqfx + qfx;
                                elseif upscb(i).flux(5) < 0
                                    qfx = -upscb(i).flux(4) * upscb(i).flux(6);
                                    tqfx = tqfx + qfx;
                                end
                                qfy = abs(upscb(i).flux(7)) * upscb(i).flux(4);
                                tqfy = tqfy + qfy;
                                %                            disp(upscb(i));
                                %                            disp(tqfx);
                            end
                        end
                    end
                end
                fclose(fid5);
                fclose(fid7);
                
                qx = tqfx + tqmx;
                qy = tqfy + tqmy;
                
                kxy = (qx * 0.001) / dx;
                kyy = (qy * 0.001) / dx;
                
                %            Deallocate memory
                clear upsca upscb;
                qx = 0;
                qy = 0;
                
                
                % Write results to file
                % Write results to file
                if kxx == 0 || kxx < 0 || kxx > 1
                    kxx = matrixp;
                end
                
                if kyy == 0 || kyy < 0 || kyy > 1
                    kyy = matrixp;
                end
                
                if kxy < -1 || kxy > 1 || abs(kxy)<matrixp
                    kxy = 0;
                end
                
                if kyx < -1 || kyx > 1 || abs(kyx)<matrixp
                    kyx = 0;
                end
                
                fprintf(fid20000, '%s %0.8e %0.8e %0.8e %0.8e\n', frac, kxx, kxy, kyx, kyy);
                
            end
            
        end
        
        
    else
        fprintf(fid20000, '%s %0.8e %0.8e %0.8e %0.8e\n', cf, kxx, kxy, kyx, kyy);
        fprintf(fid10000, 'grid %d does not exist.\n', num_grid);
        
        ia = 0;
        ib = 0;
        ic = 0;
        id = 0;
        fa = 0;
        fb = 0;
        tqmx = 0;
        qmx = 0;
        tqmy = 0;
        qmy = 0;
        tqfx = 0;
        qfx = 0;
        tqfy = 0;
        qfy = 0;
        qx = 0;
        qy = 0;
        kxx = 0;
        kyx = 0;
        kxy = 0;
        kyy = 0;
        matrixp = 0;
        
    end
end

fclose(fid10000);
fclose(fid20000);










