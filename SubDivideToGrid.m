% Part of package: UpsFrac
% Author: Tao Chen et al.
% Email: chentao9330@gmail.com
% Copyright (c) 2024 UpsFrac Developers.
% All rights reserved.
% Updated: Oct 2024

clear all;
clc;

%%
% Set the values for lx, ly, dx, dy, and tfrac
% Define the parameters for the model
lx = 1000;        % Length of the model in the x-direction (m)
ly = 1000;        % Length of the model in the y-direction (m)
dx = 100;         % Length of the coarse grid in the x-direction (m)
dy = 100;         % Length of the coarse grid in the y-direction (m)
tfrac = 203;      % Total number of fractures in the model
move_ratio = 0.001; % Precision factor for adjusting the corner points of fractures to enhance meshing efficiency
num_real = 1;     % Number of DFM realizations


%%
files = dir('resu_0*');
for i = 1:length(files)
    delete(files(i).name);
end

% Define the 'subr' structure in MATLAB
subr = struct('nsub', 0, 'cor', zeros(2, 4), 'a', []);

% Define other variables
i = 0; j = 0; k = 0; l = 0; m = 0; nx = 0; ny = 0;
tsub = 0; n = 0; flag = 0; maxx = 0; maxy = 0; pd = 0; plix = 0; pliy = 0;
plixy = 0; cmidp = 0; NumPairs = 0; LastSwap = 0; srar = 0; srac = 0;
sppox = 0; sppoy = 0; p = 0;
px = 0; py = 0; q = 0; mid = 0;  xcor = 0; ycor = 0; xcorc = 0; ycorc = 0; spc = 0; spd = 0;

% Allocate memory for arrays
lix = zeros(2, 2); % Example initialization for lix
liy = zeros(2, 2); % Example initialization for liy
lixy = zeros(2, 2); % Example initialization for lixy
midp = zeros(2, 2); % Example initialization for midp
sp = zeros(2, 2); % Example initialization for sp
cse = zeros(1, 1); % Example initialization for cse
csw = zeros(1, 1); % Example initialization for csw
cne = zeros(1, 1); % Example initialization for cne
cnw = zeros(1, 1); % Example initialization for cnw

% Allocate memory for 'sr' structure array
num_sr = 10; % Example number of 'subr' structures
sr(num_sr) = subr;


for reali = 1:num_real
    
    % Convert reali to a string with 2 digits
    seq = sprintf('%d', reali);
    
    % Open the file for reading
    fileID = fopen(['frac_0_' seq '.txt'], 'r');
    
    % Initialize variables
    srar = 4;
    sp = zeros(srar + 1, tfrac);
    
    % Read data from the file and perform manipulations
    for i = 1:tfrac
        data = fscanf(fileID, '%f', 5);
        sp(:, i) = data;
        
        % When k=+max, make y1 less than y2
        if sp(3, i) - sp(1, i) == 0
            if sp(2, i) > sp(4, i)
                spc = sp(4, i);
                sp(4, i) = sp(2, i);
                sp(2, i) = spc;
            end
        end
        
        % When x1 > x2, change position
        if sp(1, i) > sp(3, i)
            spc = sp(3, i);
            sp(3, i) = sp(1, i);
            sp(1, i) = spc;
            
            spd = sp(4, i);
            sp(4, i) = sp(2, i);
            sp(2, i) = spd;
        end
    end
    
    % Close the file
    fclose(fileID);
    
    % Display the 'sp' array
    % disp(sp);
    
    % Calculate nx, ny, and tsub
    nx = lx / dx;
    ny = ly / dy;
    tsub = nx * ny;
    
    % Allocate memory for arrays
    cse = zeros(1, tsub);
    csw = zeros(1, tsub);
    cne = zeros(1, tsub);
    cnw = zeros(1, tsub);
    sr = repmat(struct('a', zeros(srar + 1, tfrac)), 1, tsub);
    
    % Allocate memory for 'a' in 'sr'
    for i = 1:tsub
        sr(i).a = zeros(srar + 1, tfrac);
    end
    
    
    for j=1:tfrac
        
        maxx = nx + 1;
        maxy = ny + 1;
        
        for m = 1:ny
            for l = 1:nx
                n = l + (m-1) * nx;
                sr(n).nsub = n;
                sr(n).cor(1,1) = (l-1) * dx;
                sr(n).cor(2,1) = (m-1) * dy;
                sr(n).cor(1,2) = l * dx;
                sr(n).cor(2,2) = (m-1) * dy;
                sr(n).cor(1,3) = (l-1) * dx;
                sr(n).cor(2,3) = m * dy;
                sr(n).cor(1,4) = l * dx;
                sr(n).cor(2,4) = m * dy;
            end
        end
        
        
        if sp(3,j) < sp(1,j)
            sppox = sp(3,j);
            sp(3,j) = sp(1,j);
            sp(1,j) = sppox;
            sppoy = sp(4,j);
            sp(4,j) = sp(2,j);
            sp(2,j) = sppoy;
        end
        
        ssp(1,1) = sp(1,j);
        ssp(2,1) = sp(2,j);
        
        
        flag = 1;
        for i = 1:tsub
            csw(i) = sr(i).cor(1,1);
            cse(i) = sr(i).cor(1,4);
            cnw(i) = sr(i).cor(2,1);
            cne(i) = sr(i).cor(2,4);
            
            if ssp(1,1) >= csw(i) && ssp(1,1) <= cse(i) && ssp(2,1) >= cnw(i) && ssp(2,1) <= cne(i)
                ssp(1,1) = csw(i);
                ssp(2,1) = cnw(i);
            end
        end
        
        kl = (sp(4,j) - sp(2,j)) / (sp(3,j) - sp(1,j));
        
        xcor = 0;
        for i = 1:maxx
            xcorc = ssp(1,1) + i*dx - 0.0001;
            if xcorc < sp(3,j)
                xcor = xcor + 1;
            else
                break;
            end
        end
        
        if kl > 0
            ycor = 0;
            for i = 1:maxy
                ycorc = ssp(2,1) + i*dy;
                if ycorc < sp(4,j)
                    ycor = ycor + 1;
                else
                    break;
                end
            end
        elseif kl < 0
            ycor = 0;
            for i = 1:maxy
                ycorc = ssp(2,1) + dy - i*dy;
                if ycorc > sp(4,j)
                    ycor = ycor + 1;
                else
                    break;
                end
            end
        end
        
        pd = 2;
        plix = 2 + xcor;
        pliy = ycor;
        plixy = 2 + xcor + ycor;
        
        lix = zeros(pd, plix);
        liy = zeros(pd, pliy);
        lixy = zeros(pd, plixy);
        cmidp = plixy - 1;
        midp = zeros(pd, cmidp);
        
        lix = zeros(2, 2);
        lix(1, 1) = sp(1, j);
        lix(2, 1) = sp(2, j);
        lix(1, 2) = sp(3, j);
        lix(2, 2) = sp(4, j);
        
        if (sp(3, j) - sp(1, j) == 0)
            for i = 1:pliy
                py = ssp(2, 1) + i * dy;
                liy(2, i) = py;
                liy(1, i) = sp(1, j);
            end
            
            for k = 1:plixy
                if k <= plix
                    lixy(1, k) = lix(1, k);
                    lixy(2, k) = lix(2, k);
                else
                    lixy(1, k) = liy(1, k - plix);
                    lixy(2, k) = liy(2, k - plix);
                end
            end
            
            NumPairs = plixy - 1;
            
            while true
                if NumPairs == 0
                    break;
                end
                
                LastSwap = 1;
                for i = 1:NumPairs
                    if lixy(2, i) > lixy(2, i + 1)
                        Tempx = lixy(2, i);
                        lixy(2, i) = lixy(2, i + 1);
                        lixy(2, i + 1) = Tempx;
                        
                        Tempy = lixy(1, i);
                        lixy(1, i) = lixy(1, i + 1);
                        lixy(1, i + 1) = Tempy;
                        
                        % Record position of last swap
                        LastSwap = i;
                    end
                end
                
                NumPairs = LastSwap - 1;
            end
            
        else
            if plix > 2
                for i = 1:plix-2
                    px = ssp(1, 1) + i * dx;
                    lix(1, i + 2) = px;
                    lix(2, i + 2) = kl * (lix(1, i + 2) - sp(1, j)) + sp(2, j);
                end
            end
            
            if kl > 0
                for i = 1:pliy
                    py = ssp(2, 1) + i * dy;
                    liy(2, i) = py;
                    liy(1, i) = (liy(2, i) - sp(2, j)) / kl + sp(1, j);
                end
            elseif kl < 0
                for i = 1:pliy
                    py = ssp(2, 1) + dy - i * dy;
                    liy(2, i) = py;
                    liy(1, i) = (liy(2, i) - sp(2, j)) / kl + sp(1, j);
                end
            end
            
            for k = 1:plixy
                if k <= plix
                    lixy(1, k) = lix(1, k);
                    lixy(2, k) = lix(2, k);
                else
                    lixy(1, k) = liy(1, k - plix);
                    lixy(2, k) = liy(2, k - plix);
                end
            end
            
            NumPairs = plixy - 1;
            
            while true
                if NumPairs == 0
                    break;
                end
                
                LastSwap = 1;
                for i = 1:NumPairs
                    if lixy(1, i) > lixy(1, i + 1)
                        Tempx = lixy(1, i);
                        lixy(1, i) = lixy(1, i + 1);
                        lixy(1, i + 1) = Tempx;
                        
                        Tempy = lixy(2, i);
                        lixy(2, i) = lixy(2, i + 1);
                        lixy(2, i + 1) = Tempy;
                        
                        % Record position of last swap
                        LastSwap = i;
                    end
                end
                
                NumPairs = LastSwap - 1;
            end
        end
        
        
        % Initialize variables
        midp = zeros(2, cmidp);
        sr_a = zeros(4, tsub);
        
        % Calculate midpoints
        for i = 1:cmidp
            midp(1, i) = (lixy(1, i) + lixy(1, i+1)) * 0.5;
            midp(2, i) = (lixy(2, i) + lixy(2, i+1)) * 0.5;
        end
        
        %         % Process data
        for p = 1:cmidp
            for i = 1:tsub
                if (midp(1, p) >= csw(i) && midp(1, p) <= cse(i) && midp(2, p) >= cnw(i) && midp(2, p) <= cne(i))
                    sr(i).a(1,j) = lixy(1,p) - sr(i).cor(1,1);
                    sr(i).a(2,j) = lixy(2,p) - sr(i).cor(2,1);
                    sr(i).a(3,j) = lixy(1,p+1) - sr(i).cor(1,1);
                    sr(i).a(4,j) = lixy(2,p+1) - sr(i).cor(2,1);
                else
                    sr(i).a = zeros(size(sr(i).a));
                end
                %
                %                 % Write to file
                %                  fid = fopen(['null_', strtrim(seq), '.txt'], 'w');
                
                if sr(i).a(1,j) == sr(i).a(3,j) && sr(i).a(2,j) == sr(i).a(4,j)
                    %                      fprintf(fid, '%d\n', sr(i).nsub);
                    %                      fprintf(fid, '%f\n', sr(i).a(1,j));
                    %                      fprintf(fid, '%f\n', sr(i).a(2,j));
                    %                      fprintf(fid, '%f\n', sr(i).a(3,j));
                    %                      fprintf(fid, '%f\n', sr(i).a(4,j));
                    
                elseif sqrt((sr(i).a(1,j) - sr(i).a(3,j))^2 + (sr(i).a(2,j) - sr(i).a(4,j))^2) < 0.1
                    % delete short fractures
                    
                else
                    if sr(i).a(1,j) == 0 && sr(i).a(2,j) == 0
                        sr(i).a(1,j) = dx/1000 + j * move_ratio;
                        
                    elseif sr(i).a(1,j) == 0 && (dy - sr(i).a(2,j)) == 0
                        sr(i).a(1,j) = dx/1000 + j * move_ratio;
                        
                    elseif (2 - sr(i).a(1,j)) == 0 && sr(i).a(2,j) == 0
                        sr(i).a(1,j) = dx-dx/1000 - j * move_ratio;
                        
                    elseif (2 - sr(i).a(1,j)) == 0 && (dy - sr(i).a(2,j)) == 0
                        sr(i).a(1,j) = dx-dx/1000 - j * move_ratio;
                        
                    elseif sr(i).a(3,j) == 0 && sr(i).a(4,j) == 0
                        sr(i).a(3,j) = dx/1000 + j * move_ratio;
                        
                    elseif sr(i).a(3,j) == 0 && (dy - sr(i).a(4,j)) == 0
                        sr(i).a(3,j) = dx/1000 + j * move_ratio;
                        
                    elseif (dx - sr(i).a(3,j)) == 0 && sr(i).a(4,j) == 0
                        sr(i).a(3,j) = dx-dx/1000 - j * move_ratio;
                        
                    elseif (dx - sr(i).a(3,j)) == 0 && (dy - sr(i).a(4,j)) == 0
                        sr(i).a(3,j) = dx-dx/1000 - j * move_ratio;
                        
                    end
                    
                    fid_output = fopen(['resu_0_', strtrim(seq), '.txt'], 'a');
                    fprintf(fid_output, '%d %f %f %f %f %e\n', sr(i).nsub, sr(i).a(1,j), sr(i).a(2,j), sr(i).a(3,j), sr(i).a(4,j), sp(5,j));
                    % fprintf(fid_output, '%d\n', j); % Uncomment if needed
                    
                    fclose(fid_output);
                end
                
                %                  fclose(fid);
                
            end
        end
    end
end



