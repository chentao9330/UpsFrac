% Part of package: UpsFrac
% Author: Tao Chen et al.
% Email: chentao9330@gmail.com
% Copyright (c) 2024 UpsFrac Developers.
% All rights reserved.
% Updated: Oct 2024



clear; clc;

elines = load('epcomf.txt');
kxx = elines(:, 2);
kxy = elines(:, 3);
kyx = elines(:, 4);
kyy = elines(:, 5);

% 对称张量
skxy = 0.5 * (kxy + kyx);

% �? kxx 重塑�? 10x10 网格
kxx_grid = reshape(kxx, [10, 10])';
skxy_grid = reshape(skxy, [10, 10])';
kyy_grid = reshape(kyy, [10, 10])';

% 创建图形窗口
figure('Units', 'inches', 'Position', [1, 1, 6, 0.94]);

% 第一个子�?
subplot(1, 3, 1);
kxx_grid = log10(kxx_grid);
imagesc(kxx_grid);
colorbar; % 添加颜色�?caxis([1e-12,1e-16]); 
xlabel('x');
ylabel('y');
title('kxx with Logarithmic Color Scale');
set(gca, 'YDir', 'normal'); % 纵轴从下�?上增�?
axis equal;
set(gca, 'XTick', 1:10);


% 第二个子�?
subplot(1, 3, 2);
skxy_grid = log10(abs(skxy_grid));
imagesc(skxy_grid);
colorbar; % 添加颜色�?
xlabel('x');
ylabel('y');
title('kxy with Logarithmic Color Scale');
set(gca, 'YDir', 'normal'); % 纵轴从下�?上增�?
axis equal;
set(gca, 'XTick', 1:10);

% 第三个子�?
subplot(1, 3, 3);
kyy_grid = log10(kyy_grid);
imagesc(kyy_grid);
colorbar; % 添加颜色�?
xlabel('x');
ylabel('y');
title('kyy with Logarithmic Color Scale');
set(gca, 'YDir', 'normal'); % 纵轴从下�?上增�?
axis equal;
set(gca, 'XTick', 1:10);
