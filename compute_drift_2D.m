function [shifteddata, drift_info] = compute_drift_2D(data, specs, timings)

% COMPUTE_DRIFT_2D. Corrects for the drift of a 2D localization microscopy
% dataset. Input data are a list of localization coordinates in the format
% [x y frame]. The data are divided into temporal bins according to the 
% inputted frame values. The function calculates the displacement between
% bins using a mean shift algorithm. It relies on "Fcrosspairs.c", which
% must be compiled into a MEX file on the user's system.
%
% INPUTS:
%   DATA: 
%       a list of localization coordinates in the format [x y frame]
%   SPECS:
%       parameters and options for the correction. Fields can be set by
%       calling drift_default_meanshift for default values and adjusting
%       accordingly.
%   TIMINGS:
%       real times for each frame. An optional argument.
%    
% OUTPUTS:
%   SHIFTEDDATA: 
%       the drift corrected dataset
%   DRIFT_INFO: 
%       includes the drift curve, interpolated values for the drift, and
%       error estimates for each pairwise shift calculation  


% Input: data: A nx3 matrix with columns x, y, and frame
% specs: options for the alignment. Fields:
% npoints_for_alignment. The number of temporal bins into which the data
% are placed.
% nframes_per_alignment. The number of frames of data used for each alignment.
% rmax. The maximum pairwise distance considered.
% outlier_error. The maximum acceptable error.
% Output: shifteddata. The drift corrected-data.
% drift_info. Includes the drift curves at each time bin and other
% information about the alignment.

% sort by frame
olddata = data;
data = sort(data, 3);
nframes = max(data(:, 3));

% Get parameters
nTimeBin = round(specs.nTimeBin);
binwidth = floor(nframes/nTimeBin);

rmax = specs.rmax;
delta_broad = specs.delta_broad;
calc_error = specs.calc_error;
broadsweep = specs.broadsweep;
interp_method = specs.interp_method;
outlier_error = specs.outlier_error;
maxiter = specs.maxiter;

roi_width = max(data(:,1)) - min(data(:,1));
roi_height = max(data(:,2)) - min(data(:,2));
area = roi_width * roi_height;

% Compute which frames belong to which time bins
binspacing = (nframes - binwidth)/(nTimeBin-1);
firstframe = round(1 + (0:nTimeBin-1)*binspacing);
lastframe = min(nframes, firstframe + binwidth - 1);

% initialize variables
xshift = zeros(nTimeBin, nTimeBin);
yshift = zeros(nTimeBin, nTimeBin);

dxshift = zeros(nTimeBin, nTimeBin);
dyshift = zeros(nTimeBin, nTimeBin);

ntruepairs = zeros(nTimeBin, nTimeBin);
nfalsepairs = zeros(nTimeBin, nTimeBin);
loc_error = zeros(nTimeBin, nTimeBin);
iter_toconverge = zeros(nTimeBin, nTimeBin);

% helper function for getting data
getnthdata = @(n) data(firstframe(n)<=data(:,3) & data(:,3)<=lastframe(n), :);
refdata = getnthdata(1);

% correlate reference data to get a sense for how long the molecules
% stay on, as well as for the appropriate noutmax
x1 = refdata(:,1); y1 = refdata(:,2); t1 = refdata(:,3);
x2 = x1; y2 = y1; t2 = t1;
taumax = max(t1); 
noutmax = 1e8;

[dx, dy] = crosspairs(x1, y1, t1, x2, y2, t2, rmax, 1, taumax, noutmax);
noutmax = 2*length(dx); % We'll use this later. This is an upper bound for the number of pairs.

xstart = 0; ystart = 0;
[xshift_broadsweep, yshift_broadsweep] = meanshift_2D(dx, dy, xstart, ystart, delta_broad, maxiter);
rho = length(x1)^2 / area;
[~, ~, ntruepairs_broad, ~, loc_error_auto] = meanshift_2D_error(dx, dy, xshift_broadsweep(end), yshift_broadsweep(end), delta_broad, rmax, rho);
loc_prec = loc_error_auto/sqrt(2);
delta_narrow = loc_prec*specs.delta_narrow_ratio;
mean_nframeson = (ntruepairs_broad/2)/length(x1) + 1; % relevant to the error calculation, too

% can use for loop instead of parfor if no parallel capabilities are available
parfor i = 1:nTimeBin - 1
    refdata  = getnthdata(i);
    
    xshift_temp = zeros(nTimeBin, 1); 
    yshift_temp = zeros(nTimeBin, 1);
    dxshift_temp = zeros(nTimeBin, 1); 
    dyshift_temp = zeros(nTimeBin, 1);
    
    ntruepairs_temp = zeros(nTimeBin, 1); nfalsepairs_temp = zeros(nTimeBin, 1);
    loc_error_temp = zeros(nTimeBin, 1);
    iter_toconverge_temp = zeros(nTimeBin, 1);
    for j = i+1:nTimeBin
        thisdata = getnthdata(j);
        
        x1 = refdata(:,1); y1 = refdata(:,2); t1 = refdata(:,3);
        x2 = thisdata(:,1); y2 = thisdata(:,2); t2 = thisdata(:,3);
                
        % translate this data by the previous shifts
        xoffset = xshift_temp(j-1); yoffset = yshift_temp(j-1);     
        x2 = x2 - xoffset; y2 = y2 - yoffset;
        taumin = min(min(t1) - max(t2), min(t1) - max(t2)); taumax = -taumin;
 
        [dx, dy] = crosspairs(x1, y1, t1, x2, y2, t2, rmax, taumin, taumax, noutmax);
        
        xstart = 0; ystart = 0;
        if broadsweep % necessary if there are large shifts between consecutive bins
            [xshift_broadsweep, yshift_broadsweep] = meanshift_2D(dx, dy, xstart, ystart, delta_broad, maxiter);
            xstart = xshift_broadsweep(end); ystart = yshift_broadsweep(end);
            if isnan(xstart)
                xstart = 0; ystart = 0;
            end
        end
        
        [xshift_narrow, yshift_narrow] = meanshift_2D(dx, dy, xstart, ystart, delta_narrow, maxiter);
        iter_toconverge_temp(j) = length(xshift_narrow) - 1;
        xshift_temp(j) = xshift_narrow(end) + xoffset;
        yshift_temp(j) = yshift_narrow(end) + yoffset;
        
        if isnan(xshift_temp(j))
            xshift_temp(j) = xoffset;
            yshift_temp(j) = yoffset;
        end
        
        if calc_error
            rho = length(x1) * length(x2) / area;
            [dxshift_temp(j), dyshift_temp(j), ntruepairs_temp(j), nfalsepairs_temp(j), loc_error_temp(j)] = meanshift_2D_error(dx, dy, ...
                xshift_narrow(end), yshift_narrow(end), delta_narrow, rmax, rho, mean_nframeson, loc_prec);
            if isnan(dxshift_temp(j))
                dxshift_temp(j) = rmax;
                dyshift_temp(j) = rmax;
            end
        end
    end
    
    % updating the matrices with the temporary rows
    xshift(i, :) = xshift_temp; yshift(i, :) = yshift_temp;
    dxshift(i, :) = dxshift_temp; dyshift(i, :) = dyshift_temp;
    ntruepairs(i, :) = ntruepairs_temp; nfalsepairs(i, :) = nfalsepairs_temp;
    loc_error(i, :) = loc_error_temp;
    iter_toconverge(i, :) = iter_toconverge_temp;
%     fprintf('Finished with time bin %i\n', i)
end

% check if the median of the predicted error is > loc_prec/4. If so, some of the
% alignments may be false peaks

med_error = median(dxshift(dxshift > 0));
if med_error > loc_prec/4
    warning('Errors between alignments are large: consider decreasing the number of time bins.')
end

% redundant calculation is adapted from supplementary software in "Wang et al."

% total number of iterations is nTimeBin*(nTimeBin-1)/2
nelements = (nTimeBin^2-nTimeBin)/2;
A = zeros(nelements, 2);
R = zeros(nelements, 2);
w = zeros(nelements, 1);
count = 1;
for i=1:nTimeBin-1
    for j=i+1:nTimeBin
        A(count, i) = -1;
        A(count, j) = 1;
        R(count, 1) = xshift(i,j);
        R(count, 2) = yshift(i,j);
        if calc_error
            w(count, 1) = 1/(dxshift(i,j))^2;
        else
            w(count, 1) = 1;
        end
        count = count+1;
    end
end

A = A(:, 2:end);

[D, stderr, mse] = lscov(A, R, w);
err=A*D-R;
b=R;
rowerr = zeros(size(A,1),2);
for i=1:size(A,1)
    rowerr(i,1) = sqrt(err(i,1)^2+err(i,2)^2);
end
rowerr(:,2)=1:size(A,1);
rowerr=flipud(sortrows(rowerr,1));

index=rowerr(find(rowerr(:,1)>outlier_error),2);
noutliers = numel(index);
fraction_outliers = noutliers/size(A,1);

% remove outliers while preventing rank deficiency of A
for i=1:size(index,1)
    flag = index(i);
    tmp=A;
    tmp(flag,:)=[];
    if rank(tmp,1)==(nTimeBin-1)
        A(flag,:)=[];
        b(flag,:)=[];
        w(flag,:)=[];
        sindex=find(index>flag);
        index(sindex)=index(sindex)-1;
    else
        tmp=A;
    end
end

% A(index,:)=[];
% b(index,:)=[];
% w(index,:)=[];

[D, stderr, mse] = lscov(A, b, w);

if nargin < 3
    timings = 1:nframes;
end

midtiming = zeros(nTimeBin, 1);
midtiming_frames = zeros(nTimeBin, 1);
for i = 1:nTimeBin
    midtiming(i) = mean(timings(firstframe(i):lastframe(i)));
    midtiming_frames(i) = mean(firstframe(i):lastframe(i));
end

drift = vertcat([0,0], D);
stderr = vertcat([0,0], stderr);
xinterp = interp1(midtiming, drift(:,1), timings, interp_method, 'extrap');
yinterp = interp1(midtiming, drift(:,2), timings, interp_method, 'extrap');

% info and diagnostics
drift_info.xshift = xshift;
drift_info.yshift = yshift;
drift_info.dxshift = dxshift;
drift_info.dyshift = dyshift;
drift_info.xinterp = xinterp;
drift_info.yinterp = yinterp;
drift_info.binwidth = binwidth;
drift_info.timings = timings;
drift_info.midtiming = midtiming;
drift_info.midtiming_frames = midtiming_frames;

drift_info.rowerr = rowerr;
drift_info.drift = drift;
drift_info.stderr = stderr;
drift_info.mse = mse;
drift_info.noutliers = noutliers;
drift_info.fraction_outliers = fraction_outliers;

drift_info.ntruepairs = ntruepairs;
drift_info.nfalsepairs = nfalsepairs;
drift_info.loc_error = loc_error;
drift_info.iter_toconverge = iter_toconverge;

drift_info.delta_narrow = delta_narrow;
drift_info.loc_prec = loc_prec; 
drift_info.med_error = med_error;

% Shift the data by the interpolated drift values at each time point.
shifteddata = olddata;
for i = 1:nframes
    index = find(data(:, 3) == i);
    shifteddata(index, 1) = olddata(index, 1) - xinterp(i);
    shifteddata(index, 2) = olddata(index, 2) - yinterp(i);
end
shifteddata(:,3) = olddata(:,3);