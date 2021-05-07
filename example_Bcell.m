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
driftspecs.delta_broad = 100; % radius of consideration for broad sweep
driftspecs.delta_narrow_ratio = 3; % delta is set to be this times the localization precision
driftspecs.interp_method = 'linear'; % for interpolating drift by frame
driftspecs.outlier_error = 50; % remove points with greater error than this value

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
