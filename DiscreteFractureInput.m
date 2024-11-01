% Part of package: UpsFrac
% Author: Tao Chen et al.
% Email: chentao9330@gmail.com
% Copyright (c) 2024 UpsFrac Developers.
% All rights reserved.
% Updated: Oct 2024

%%
% Set the number of DFM realizations
num_real = 1;     % Number of DFM realizations



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for reali=1:num_real
    
    seq=num2str(reali)
    realf=strcat('reali_',seq)
    
    
    mkdir (realf)
    resu_file=['resu_0_',seq,'.txt'];
    a=textread(resu_file);
    [r,c]=size(a)
    pathr=pwd
    
    fn=r
    sub=100
    cd (realf)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    path=pwd
    % creat fracture files for DFM simulation and Computation effective permeability
    
    for j=1:fn
        
        for i=1:sub
            
            if a(j,1)==i
                num=i
                t=num2str(num)
                
                %difference between windows (\) and mac (/)
                subf=strcat(path,'/',t,'.txt')
                fidd = fopen(subf,'a');
                
                mkdir RESULT
                subu=strcat('./RESULT/',t,'f.txt')
                fidu = fopen(subu,'a');
                
                ps=a(j,2:3)
                pe=a(j,4:5)
                aperture=a(j,6)
                fprintf(fidd,'%f\t%f\t;%f\t%f\t;\n',ps,pe)
                fprintf(fidu,'%f\t%f\t%f\t%f\t%g\n',ps,pe,aperture)
                %fprintf(fidu,'%f\t%f\t%f\t%f\t\%g\t\n',ps,pe,aperture)
                fclose('all')
            end
            
        end
        
    end
    
    %read standard DFM imput file
    % for windows:
    %subm=strcat(pathr,'\matlabFractureGrid_s1.m')
    % for mac:
    subm=strcat(pathr,'/InputTemplate_tpfa.m')
    
    fid = fopen(subm,'r');
    tline = fgetl(fid);
    n=1;
    header_sh={};
    while ischar(tline)
        if n<=227
            header_sh{n}=tline;
        end;
        
        tline = fgetl(fid);
        n=n+1;
    end;
    
    fclose(fid);
    
    %read creat DFM imput files for each subregion which contians fracture(s)
    
    for SubIndex=1:sub
        
        for j=1:fn
            
            if a(j,1)==SubIndex
                
                str=num2str(SubIndex)
                % for windows:
                %file_out=strcat(path,'\s',str,'.m')
                % for mac:
                file_out=strcat(path,'/s',str,'.m')
                
                fid = fopen(file_out,'wt');
                
                % copy line 1-22 from the template file 'matlabFractureGrid_s1.m'
                for i=1:22
                    fprintf(fid,'%s\n',header_sh{i})
                end
                
                % input fracture nodels based on the template file 'matlabFractureGrid_s1.m'
                for j=1:fn
                    if a(j,1)==SubIndex
                        num=SubIndex
                        ps=a(j,2:3)
                        pe=a(j,4:5)
                        fprintf(fid,'%.3f\t%.3f\t;%.3f\t%.3f\t;\n',ps,pe)
                    end
                    
                end
                
                
                fprintf(fid,'%s\n',header_sh{24})
                
                fprintf(fid,'subtxt=''%d.txt''\n',SubIndex)
                
                
                % copy line 25-227 in the template file 'matlabFractureGrid_s1.m'
                
                for i=25:56
                    fprintf(fid,'%s\n',header_sh{i})
                end
                
                mat_aper=strcat('leng_aper=dlmread(''RESULT/',str,'f.txt'')')
                fprintf(fid,'%s\n',mat_aper)
                
                for i=57:227
                    fprintf(fid,'%s\n',header_sh{i})
                end
                
                fclose(fid)
                
            end
        end
    end
    
    delete *.txt
    %movefile('RESULT', 'RESULT1');
    
    cd ..
    
end

