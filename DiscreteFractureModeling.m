% Alghalandis Discrete Fracture Network Engineering (ADFNE)
% Author: Younes Fadakar Alghalandis
% Copyright (c) 2016 Alghalandis Computing @ http://alghalandis.net
% All rights reserved.
%
% Part of package: UpsFrac
% Modified by: Tao Chen et al.
% Email: chentao9330@gmail.com
% Copyright (c) 2024  UpsFrac Developers.
% All rights reserved.
% Modification Date: Oct 2024
% Description of modifications:
%  (1) power law length distribution;
%  (2) correlating aperture to fracture length;
%  (3) combining stochastic and determined fractures;
%  (4) export DFN in files;

clear all
clc

%% Paremeters of the fracture network
%number of the fractures / density
n = 200;

%length of the square region
L = 1000;

%%
% Set the number of DFM realizations
num_real = 1;     % Number of DFM realizations

for i = 1:num_real
    a=num2str(i)
    
    %% region of the fracture network
    rgn = [0,L,0,L];
    
    %% modeling location of the fractures
    pts = L*rand(n,2);
    scatter(pts(:,1),pts(:,2))
    axis([0 L 0 L])
    
    %% modeling orientation of the fracture network
    %orientation parameter
    kappa = 0;%random when kappa=0
    theta=3*pi/4;
    ags = circ_vmrnd(theta,kappa,n);
    
    %% modeling length of the fractures
    %power-law distribution
    % Parameters
    l_min = 100;
    l_max = 1000;
    number = n;
    exponent = 2.5;
    
    % Generate fracture lengths using the power-law function
    lhs = rndm_powerlaw(l_min, l_max, exponent*(-1)+1, number);
    
    %% Generate apertures: length (l)-aperture(a) coorelation: a=A*l(exp(b)),
    % b=0, constant;b=0.5, sublinear;b=1.0, linear;
    A=2.3e-5;
    b=0.5;
    Aperture_0=A*lhs.^b;
    
    %% generation of the fracture network
    [dx,dy] = pol2cart(ags,0.5*lhs);
    olines = [pts(:,1)-dx,pts(:,2)-dy,pts(:,1)+dx,pts(:,2)+dy]; %original
    lines = ClipLines2D(olines,rgn);
    

    %% combine stochastic and determintic fractures
    determ_data = load('frac_determintic.txt');
    if ~isempty(determ_data)
        lines = [lines; determ_data(:,1:4)];
        Aperture_0 = [Aperture_0; determ_data(:,5)]; 
    end
    
    %% Visualisation
    clf
    [X,Y] = LinesToXYnan2D(lines);
    
    % Assuming X, Y, and a are defined
    numLines = size(X, 2); % Number of lines
    colors = zeros(numLines, 3); % Initialize color array
    
    % Normalize 'a' to range [0, 1]
    a_normalized = (Aperture_0 - min(Aperture_0)) / (max(Aperture_0) - min(Aperture_0));
    
    % Create a colormap from blue to red
    for i = 1:numLines
        colors(i, :) = [a_normalized(i), 0, 1 - a_normalized(i)]; % Interpolate color
    end
    
    % Plot each line with the corresponding color
    hold on; % Ensure all lines are plotted on the same figure
    for i = 1:numLines
        for j=1:4
            if lines(i,j)<0
                lines(i,j)=0
            end
            if lines(i,j)>L
                lines(i,j)=L
            end
        end
        plot(X(:, i), Y(:, i), 'Color', colors(i, :)); % Use color based on 'a'
    end
    
    axis([0 L 0 L])
    axis equal;
    
    % Set the plot range from 0 to L for both X and Y axes
    xlim([0 L]);
    ylim([0 L]);
    
    grid on;
    set(gca,'Xtick',0:100:L);
    set(gca,'Ytick',0:100:L);
    
    box on;
    
    %% save the modeling results
    
    % fid=fopen(['frac_dfm',a,'.txt'],'w');
    % fprintf(fid,'%8.6f %8.6f; %8.6f %8.6f; \n',lines')
    %
    % fid=fopen(['frac_',a,'.txt'],'w');
    % fprintf(fid,'%8.6f %8.6f %8.6f %8.6f \n',lines')
    
    lines(:,5)=Aperture_0
    
    fid=fopen(['frac_0_',a,'.txt'],'w');
    fprintf(fid,'%8.6f %8.6f %8.6f %8.6f %6.2g\n',lines')
    
end