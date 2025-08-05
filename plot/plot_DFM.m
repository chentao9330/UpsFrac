% Part of package: UpsFrac
% Author: Tao Chen et al.
% Email: chentao9330@gmail.com
% Copyright (c) 2024 UpsFrac Developers.
% All rights reserved.
% Updated: Oct 2024


L=1000;
lines = load( 'frac_0_1.txt');


figure;

[X,Y] = LinesToXYnan2D(lines(:,1:4));
Aperture_0=lines(:,5);

% Assuming X, Y, and a are defined
numLines = size(X, 2); % Number of lines
colors = zeros(numLines, 3); % Initialize color array

% Normalize 'a' to range [0, 1]
a_normalized = (Aperture_0 - min(Aperture_0)) / (max(Aperture_0) - min(Aperture_0));

% Create a colormap from blue to red
% Normalize 'a' to range [0, 1]
if (max(Aperture_0) - min(Aperture_0)) == 0
    a_normalized = 0.5 * ones(size(Aperture_0));  % constant aperture
else
    a_normalized = (Aperture_0 - min(Aperture_0)) / (max(Aperture_0) - min(Aperture_0));
end

% Plot each line with the corresponding color
hold on; % Ensure all lines are plotted on the same figure
for i = 1:numLines
        
     plot(X(:, i), Y(:, i), 'Color', colors(i, :)); % Use color based on 'a'
end

axis([0 L 0 L])
box on;
axis equal;

% Set the plot range from 0 to L for both X and Y axes
xlim([0 L]);
ylim([0 L]);


grid on;
set(gca,'Xtick',0:100:L);
set(gca,'Ytick',0:100:L);
% Set the grid color to grey