% Part of package: UpsFrac
% Author: Tao Chen et al.
% Email: chentao9330@gmail.com
% Copyright (c) 2024 UpsFrac Developers.
% All rights reserved.
% Updated: Oct 2024



%%
% Set the number of DFM realizations
num_real = 1;     % Number of DFM realizations
num_grid = 100;   % Number of Cartesian grid for a realization

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for reali=1:num_real
    
    %for reali=[1]
    seq=num2str(reali)
    realf=strcat('reali_',seq)
    
    cd (realf)
    delete('*.txt')
    
    folder  =pwd
    realmatrix=strcat('realmatrix',seq,'.txt')
    fidmatrix = fopen(realmatrix,'a');
    
    
    fileID = fopen('error_simu.txt', 'a');
    
    
    %grid number:100
    for i=1:num_grid
        grids=num2str(i)
        runf=strcat('s',grids,'.m')
        
        currentDir = pwd;
        
        runf = fullfile(currentDir, runf);
        
        %judge if the file code exist (i.e., the grid contain fractures)
        if exist(runf,'file')==2
            % DFM modeling error occurs, pass and go on simulation
            try
                newFile = runf;
                run(newFile)
                
            catch ME
                
                fprintf(fileID, 'Error occurred on grid %s\n', grids);
                fprintf('here is %s\n', grids);
            end
        else
            disp('no exist');
            fprintf(fidmatrix,'matrix permeability: %i\t;\n',i)
        end
        
    end
    
    cd ..
    close all
    
    fclose('all')
    
    cd (realf)
end

