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