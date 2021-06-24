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
% along with MEAN SHIFT DRIFT CORRECTION.  If not, see <https://www.gnu.org/licenses/>

function I_blur = render_Bcell(data)
x = data(:,1);
y = data(:,2);
z = data(:,3);
z = z + 300;

psize = 25;
xedges = min(x):psize:max(x);
yedges = min(y):psize:max(y);

[Ixy, ~, ~, binx, biny] = histcounts2(x, y, xedges, yedges); 

allz = zeros(size(Ixy));
countz = allz;
for i=1:numel(binx)
    if binx(i)>0
       % disp('.')
    allz(binx(i), biny(i)) = allz(binx(i), biny(i))+z(i);
    countz(binx(i), biny(i)) = countz(binx(i), biny(i))+1;
    end
end

blur_by = 1;
off_thresh = .05;
allz_blur = gausblur(allz,blur_by);
countz_blur = gausblur(countz,blur_by);
im2 = (allz_blur./countz_blur)'/250;
nd2 = countz_blur' < off_thresh;
rgb2 = ind2rgb(uint8(im2*255),parula(256));
di2 = ind2rgb(~nd2, [0 0 0; 1 1 1]);

I_blur = rgb2.*di2;
end