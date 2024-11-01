% Part of package: UpsFrac 
% Author: Tao Chen
% Email: chentao9330@gmail.com
% Copyright (c) 2024 Tao Chen
% All rights reserved. 
% Updated: Oct 2024
% Define the function for power-law generation
% https://stackoverflow.com/questions/31114330/python-generating-random-numbers-from-a-power-law-distribution


function xx = rndm_powerlaw(a, b, g, size)
    % Power-law gen for pdf(x)\propto x^{g-1} for a<=x<=b
    r = rand(size, 1);
    ag = a^g;
    bg = b^g;
    xx = (ag + (bg - ag) * r).^(1/g);
end