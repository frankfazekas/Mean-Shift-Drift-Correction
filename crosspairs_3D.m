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

function [dxout, dyout, dzout, dtout] = crosspairs_3D(x1, y1, z1, t1, x2, y2, z2, t2, rmax,...
        taumin, taumax, noutmax);

% Sort based on x1
[x1, s1] = sort(x1(:));
x1 = x1';
[x2, s2] = sort(x2(:));
x2 = x2';

y1 = y1(s1);
z1 = z1(s1);
t1 = t1(s1);
y2 = y2(s2);
z2 = z2(s2);
t2 = t2(s2);

[dxout, dyout, dzout, dtout] = Fcrosspairs_3D(x1, y1, z1, t1, x2, y2, z2, t2,...
                                    rmax, taumin, taumax, int64(noutmax));