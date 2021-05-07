% load in the data
data = csvread('Data/NPCdata.csv');

%% Setting options for the algorithm
% drift_default_meanshift returns reasonable default values for the options
driftspecs = drift_default_meanshift;

% modify the default options individually in the driftspecs structure
driftspecs.nTimeBin = 50; % Number of time intervals to divide the data into
driftspecs.rmax = 500; % Maximum pairwise distances extracted
driftspecs.maxiter = 200; % max iterations for the mean shift algorithm
driftspecs.calc_error = 1; % provide error estimates for mean shift
driftspecs.broadsweep = 0; % set to 1 if large shifts are expected. Not required here.
driftspecs.delta_broad = 100; % radius of consideration for 
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
