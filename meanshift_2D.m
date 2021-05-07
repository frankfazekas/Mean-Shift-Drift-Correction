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