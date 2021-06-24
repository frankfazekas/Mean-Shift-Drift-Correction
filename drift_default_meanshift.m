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

function driftspecs = drift_default_meanshift()

driftspecs.nTimeBin = 20; % number of points used for the alignment
driftspecs.rmax = 500; % maximum pairwise displacements considered
driftspecs.maxiter = 200; % max iterations for the mean shift algorithm
driftspecs.calc_error = 1; % provide error estimates for mean shift
driftspecs.broadsweep = 0; % set to 1 if large shifts are expected
driftspecs.delta_broad = 100; % radius of consideration for 
driftspecs.delta_narrow_ratio = 3; % delta is set to be this times the localization precision
driftspecs.interp_method = 'linear'; % for interpolating drift by frame
driftspecs.outlier_error = 50; % remove points with greater error than this value