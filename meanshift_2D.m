% Copyright (C) 2021 Frank Fazekas, Thomas Shaw, and Sarah Veatch
% This file is part of MEAN SHIFT DRIFT CORRECTION.
% MEAN SHIFT DRIFT CORRECTION is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% MEAN SHIFT DRIFT CORRECTION is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% You should have received a copy of the GNU General Public License
% along with MEAN SHIFT DRIFT CORRECTION.  If not, see <https://www.gnu.org/licenses/>.

function [xshift, yshift] = meanshift_2D(dx, dy, xstart, ystart, delta, maxiter)
xshift = xstart;
yshift = ystart;

n = 0;
xshift_diff = 1; % just to make loop condition work
yshift_diff = 1;

% Gaussian mean shift approach
while (xshift_diff ~= 0 && yshift_diff ~= 0 && n < maxiter && ~isnan(xshift(end)))
    % find all pairs within a distance of delta from the current shift
    r = sqrt((dx - xshift(end)).^2 + (dy - yshift(end)).^2);
    closei = find(r < delta);
    
    % Update the shift estimate
    xshift = [xshift; mean(dx(closei))];
    yshift = [yshift; mean(dy(closei))];
    
    xshift_diff = xshift(end) - xshift(end - 1);
    yshift_diff = yshift(end) - yshift(end - 1);
    
    n = n + 1; 
end