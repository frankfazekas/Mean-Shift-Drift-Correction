function [xshift, yshift, zshift] = meanshift_3D(dx, dy, dz, xstart, ystart, zstart, delta, maxiter)

xshift = xstart;
yshift = ystart;
zshift = zstart;

n = 0;
xshift_diff = 1; % just to make loop condition work
yshift_diff = 1;
zshift_diff = 1;

% Gaussian mean shift approach
while (xshift_diff ~= 0 && yshift_diff ~= 0 && zshift_diff ~= 0 && n < maxiter && ~isnan(xshift(end)))
    % there is more error in z, so we will give z displacements a discount
    % factor
    z_discount = sqrt(2);
    r = sqrt((dx - xshift(end)).^2 + (dy - yshift(end)).^2 + (dz - zshift(end)).^2/z_discount^2); % more error in z
    closei = find(r < delta);
    xshift = [xshift; mean(dx(closei))];
    yshift = [yshift; mean(dy(closei))];
    zshift = [zshift; mean(dz(closei))];
    xshift_diff = xshift(end) - xshift(end - 1);
    yshift_diff = yshift(end) - yshift(end - 1);
    zshift_diff = zshift(end) - zshift(end - 1);
    
    n = n + 1;
end