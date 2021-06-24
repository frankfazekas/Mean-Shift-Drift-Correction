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

% load in the data
data = csvread('Data/NPCdata.csv');

%% Setting options for the algorithm
% drift_default_meanshift returns reasonable default values for the options
driftspecs = drift_default_meanshift;

% modify the default options individually in the driftspecs structure
driftspecs.nTimeBin = 100; % Number of time intervals to divide the data into
driftspecs.rmax = 500; % Maximum pairwise distances extracted
driftspecs.maxiter = 200; % max iterations for the mean shift algorithm
driftspecs.calc_error = 1; % provide error estimates for mean shift
driftspecs.broadsweep = 0; % set to 1 if large shifts are expected. Not required here.
driftspecs.delta_broad = 100; % radius of consideration for the broad sweep
driftspecs.delta_narrow_ratio = 3; % delta is set to be this times the localization precision
driftspecs.interp_method = 'linear'; % for interpolating drift by frame
driftspecs.outlier_error = 50; % remove points with greater error than this value

%% Run the drift estimation and correction
[shifteddata, drift_info] = compute_drift_2D(data, driftspecs);

%% Show results
% drift_info contains the drift trajectory and various diagnostic information
% show the drift trajectory
figure; errorbar(drift_info.drift(:,1), drift_info.drift(:,2), drift_info.stderr(:,2), drift_info.stderr(:,2), drift_info.stderr(:,1), drift_info.stderr(:,1),'b-o')
xlabel('X-Drift (nm)')
ylabel('Y-Drift (nm)')
axis equal

figure; plot(drift_info.timings, drift_info.xinterp, 'b')
hold on; plot(drift_info.timings, drift_info.yinterp, 'r')
legend('X-Drift', 'Y-Drift', 'location', 'best')
xlabel('Timing')
ylabel('Drift (nm)')
axis tight

psize = 2.5;
% render the uncorrected image
x1 = data(:,1); y1 = data(:,2);
x1edges = min(x1):psize:max(x1);
y1edges = min(y1):psize:max(y1);
I1 = histcounts2(y1, x1, y1edges, x1edges);
I1 = gausblur(I1, 2*psize);
figure; imshow(3*I1)

% render the corrected image
x2 = shifteddata(:,1); y2 = shifteddata(:,2);
x2edges = min(x2):psize:max(x2);
y2edges = min(y2):psize:max(y2);
I2 = histcounts2(y2, x2, y2edges, x2edges);
I2 = gausblur(I2, 2*psize);
figure; imshow(3*I2)

%% calculate the FRC
% mask just the lower nucleus for the FRC calculation
maskx = [34072 32489 32825 37526 41268 42659 43715 41652 38582 34360 34072];
masky = [38703 41390 44556 46331 45036 43597 40430 38607 37936 38607 38703];

tokeep = zeros(size(shifteddata, 1), 1);
for i = 1:numel(tokeep)
    if inpolygon(shifteddata(i,1), shifteddata(i,1), maskx, masky)
        tokeep(i) = 1;
    end
end

shifteddata_masked = shifteddata(tokeep == 1, :);

addpath('FIREfunctions')
pixelsize = 100;
x2_pixel = shifteddata_masked(:,1)/100; y2_pixel = shifteddata_masked(:,2)/100;
x2_pixel = x2_pixel - 300; y2_pixel = y2_pixel - 300;
coords = [x2_pixel y2_pixel shifteddata_masked(:,3)];

nblocks = 31;
nreps = 1; % taking more reps slows down the computation but may improve precision
superzoom = 20;
szx = superzoom * 200;

fprintf('\n -- computing FIRE --\n')
[fire_value_px, frc_curve, fireH, fireL] = postoresolution(coords, szx, superzoom, nblocks, [], nreps);
frcplot(fire_value_px, fireL, fireH, pixelsize, superzoom, frc_curve)
fire_value = fire_value_px*pixelsize/superzoom; 
fire_value_uncertainty = (fireL-fireH)/2*pixelsize/superzoom;
