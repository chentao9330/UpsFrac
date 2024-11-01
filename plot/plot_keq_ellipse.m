% Part of package: UpsFrac
% Author: Tao Chen et al.
% Email: chentao9330@gmail.com
% Copyright (c) 2024 UpsFrac Developers.
% All rights reserved.
% Updated: Oct 2024



clear; clc;

Nx=10;
Ny=10;

elines = load( 'epcomf.txt');
kxx=elines(:,2);
kxy=elines(:,3);
kyx=elines(:,4);
kyy=elines(:,5);

% symmetrric tensor
skxy=0.5*(kxy+kyx);

% from bottom to upper plot
bxx = reshape(kxx, Nx, Ny)';
txx = bxx(end:-1:1, :); 
bkxx = reshape(txx', Nx*Ny, 1);

bxy = reshape(skxy, Nx, Ny)';
txy = bxy(end:-1:1, :); 
bkxy = reshape(txy', Nx*Ny, 1);

byy = reshape(kyy, Nx, Ny)';
tyy = byy(end:-1:1, :); 
bkyy = reshape(tyy', Nx*Ny, 1);



% Loop through each row of the matrix
for i = 1:100%Nx*Ny
  
    
    % Construct the permeability tensor
    K = [bkxx(i), bkxy(i); bkxy(i), bkyy(i)];

    
    % Calculate the eigenvalues and eigenvectors
    [eigVec, eigVal] = eig(K);
    
    % Create a grid for the ellipse
    theta = linspace(0, 2 * pi, 100);
    ellipse = [cos(theta); sin(theta)];
    
    % Scale the ellipse based on the eigenvalues
    scaledEllipse = eigVec * sqrt(eigVal) * ellipse;

    % Determine subplot position
    subplot(10, 10, i); % 10 rows, 10 columns, i-th subplot
    plot(scaledEllipse(1, :), scaledEllipse(2, :), 'b-', 'LineWidth', 2);
    axis equal; 
    
    
    maxValue = max(max(scaledEllipse));
    maxValue = abs(real(maxValue))

    xlim([-1*maxValue maxValue]);
    ylim([-1*maxValue maxValue]);
    
    

    set(gca, 'XTick', [], 'YTick', []);
    grid off;


xlabel('');
ylabel('');
title('');


box on;

hold off;

clear scaledEllipse;

end

set(gcf, 'Position', [1, 1, 800, 800]); 
subplots = findobj(gcf, 'Type', 'Axes'); 
for i = 1:length(subplots)
    set(subplots(i), 'Position', get(subplots(i), 'Position') .* [0.76, 0.72,1,1]); 
end
