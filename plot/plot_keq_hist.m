% Part of package: UpsFrac
% Author: Tao Chen et al.
% Email: chentao9330@gmail.com
% Copyright (c) 2024 UpsFrac Developers.
% All rights reserved.
% Updated: Oct 2024

%L=10;
clear; clc;

elines = load( 'epcomf.txt');
kxx=elines(:,2);
kxy=elines(:,3);
kyx=elines(:,4);
kyy=elines(:,5);

% symmetrric tensor
skxy=0.5*(kxy+kyx);



% Reshape kxx into a 10x10 grid
kxx_grid = reshape(kxx, [10, 10])';
skxy_grid = reshape(skxy, [10, 10])';
kyy_grid = reshape(kyy, [10, 10])';


% Create a figure
figure('Units', 'inches', 'Position', [1, 1, 6, 2]); 

subplot(1, 3, 1); % 1 row, 3 columns, 1st subplot


% Plot the histogram
histogram(log10(kxx));

% Add labels and title
xlabel('Value');
ylabel('Frequency');
title('Histogram of log10(kxx)');

% Optionally, set grid
grid on;

subplot(1, 3, 2); % 1 row, 3 columns, 1st subplot

% Plot the histogram
histogram(log10(abs(skxy)));

% Add labels and title
xlabel('Value');
ylabel('Frequency');
title('Histogram of log10(abs(kxy))');

% Optionally, set grid
grid on;

subplot(1, 3, 3); % 1 row, 3 columns, 1st subplot


% Plot the histogram
histogram(log10(kyy));

% Add labels and title
xlabel('Value');
ylabel('Frequency');
title('Histogram of log10(kyy)');

% Optionally, set grid
grid on;





