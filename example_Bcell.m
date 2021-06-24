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
data = csvread('Data/Bcelldata.csv');

%% Setting options for the algorithm
% drift_default_meanshift returns reasonable default values for the options
driftspecs = drift_default_meanshift;

% modify the default options individually in the driftspecs structure
driftspecs.nTimeBin = 50; % Number of time intervals to divide the data into
driftspecs.rmax = 500; % Maximum pairwise distances extracted
driftspecs.maxiter = 200; % max iterations for the mean shift algorithm
driftspecs.calc_error = 1; % provide error estimates for mean shift
driftspecs.broadsweep = 0; % set to 1 if large shifts are expected. Not required here.
driftspecs.delta_broad = 100; % radius of consideration for the broad sweep
driftspecs.delta_narrow_ratio = 3; % delta is set to be this times the localization precision
driftspecs.interp_method = 'linear'; % for interpolating drift by frame
driftspecs.outlier_error = 100; % remove points with greater error than this value

%% Run the drift estimation and correction
[shifteddata, drift_info] = compute_drift_3D(data, driftspecs);

%% Show results
% drift_info contains the drift trajectory and various diagnostic information
% show the drift trajectory
figure; errorbar(drift_info.drift(:,1), drift_info.drift(:,2), drift_info.stderr(:,2), drift_info.stderr(:,2), drift_info.stderr(:,1), drift_info.stderr(:,1),'b-o')
xlabel('X-Drift (nm)')
ylabel('Y-Drift (nm)')
axis equal

figure; errorbar(drift_info.drift(:,1), drift_info.drift(:,3), drift_info.stderr(:,3), drift_info.stderr(:,3), drift_info.stderr(:,1), drift_info.stderr(:,1),'b-o')
xlabel('X-Drift (nm)')
ylabel('Z-Drift (nm)')
axis equal

figure; plot(drift_info.timings, drift_info.xinterp, 'b')
hold on; plot(drift_info.timings, drift_info.yinterp, 'r')
hold on; plot(drift_info.timings, drift_info.zinterp, 'k')
legend('X-Drift', 'Y-Drift', 'Z-Drift', 'location', 'best')
xlabel('Timing')
ylabel('Drift (nm)')
axis tight

% render the uncorrected image
I1 = render_Bcell(data);
figure; imshow(I1)
c1 = colorbar('Ticks', [0 .4 .8 ], 'TickLabels', {'0nm', '100nm', '200nm'});
c1.Label.String = 'z (nm)';

% render the corrected image
I2 = render_Bcell(shifteddata);
figure; imshow(I2)
c2 = colorbar('Ticks', [0 .4 .8 ], 'TickLabels', {'0nm', '100nm', '200nm'});
c2.Label.String = 'z (nm)';

%% calculate the FRC
% may take a few minutes for this dataset
addpath('FIREfunctions')
x2 = shifteddata(:,1); y2 = shifteddata(:,2); z2 = shifteddata(:,3);
z2 = z2 + 1500;
pixelsize = 100;
x2_pixel = x2/100; y2_pixel = y2/100; z2_pixel = z2/100;
coords_xy = [x2_pixel y2_pixel shifteddata(:,4)];
coords_xz = [x2_pixel z2_pixel shifteddata(:,4)];

nblocks = 32;
nreps = 1; % taking more reps slows down the computation but may improve precision
superzoom = 10;
szx = superzoom * 700;

fprintf('\n -- computing FIRE --\n')
% xy FRC
[fire_value_xy_px, frc_curve_xy, fireH_xy, fireL_xy] = postoresolution(coords_xy, szx, superzoom, nblocks, [], nreps); % in super-resolution pixels
frcplot(fire_value_xy_px, fireL_xy, fireH_xy, pixelsize, superzoom, frc_curve)
fire_value_xy = fire_value_xy_px*pixelsize/superzoom; 
fire_value_xy_uncertainty = (fireL_xy-fireH_xy)/2*pixelsize/superzoom;

fprintf('\n -- computing FIRE --\n')
% xz FRC
[fire_value_xz_px, frc_curve_xz, fireH_xz, fireL_xz] = postoresolution(coords_xz, szx, superzoom, nblocks, [], nreps); % in super-resolution pixels
frcplot(fire_value_xz_px, fireL_xz, fireH_xz, pixelsize, superzoom, frc_curve)
fire_value_xz = fire_value_xz_px*pixelsize/superzoom; 
fire_value_xz_uncertainty = (fireL_xz-fireH_xz)/2*pixelsize/superzoom;