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

% Using code from "RCC.m" by Yina Wang @ Hust 2013.09.09

function [shifteddata, drift_info] = compute_drift_3D(data, specs, timings)

% COMPUTE_DRIFT_3D. Corrects for the drift of a 3D localization microscopy
% dataset. Input data are a list of localization coordinates in the format
% [x y z frame]. The data are divided into temporal bins according to the 
% inputted frame values. The function calculates the displacement between
% bins using a mean shift algorithm. It relies on "Fcrosspairs_3D.c", which
% must be compiled into a MEX file on the user's system.
%
% INPUTS:
%   DATA: 
%       a list of localization coordinates in the format [x y z frame]
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

% sort by frame
olddata = data;
data = sort(data, 4);
nframes = max(data(:, 4));

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

xshift = zeros(nTimeBin, nTimeBin);
yshift = zeros(nTimeBin, nTimeBin);
zshift = zeros(nTimeBin, nTimeBin);

dxshift = zeros(nTimeBin, nTimeBin);
dyshift = zeros(nTimeBin, nTimeBin);
dzshift = zeros(nTimeBin, nTimeBin);

ntruepairs = zeros(nTimeBin, nTimeBin);
nfalsepairs = zeros(nTimeBin, nTimeBin);
ntruepairs_z = zeros(nTimeBin, nTimeBin);
nfalsepairs_z = zeros(nTimeBin, nTimeBin);
loc_error = zeros(nTimeBin, nTimeBin);
loc_error_z = zeros(nTimeBin, nTimeBin);
iter_toconverge = zeros(nTimeBin, nTimeBin);

% helper function for getting data
getnthdata = @(n) data(firstframe(n)<=data(:,4) & data(:,4)<=lastframe(n), :);
refdata = getnthdata(1);

% correlate reference data to get a sense for how long the molecules
% stay on, as well as for the appropriate noutmax
x1 = refdata(:,1); y1 = refdata(:,2); z1 = refdata(:,3); t1 = refdata(:,4);
x2 = x1; y2 = y1; z2 = z1; t2 = t1;
taumax = max(t1); 
noutmax = 1e8;

[dx, dy, dz] = crosspairs_3D(x1, y1, z1, t1, x2, y2, z2, t2, rmax, 1, taumax, noutmax);
noutmax = 2*length(dx); % We'll use this later. This is an upper bound for the number of pairs.

xstart = 0; ystart = 0; zstart = 0;
[xshift_broadsweep, yshift_broadsweep, zshift_broadsweep] = meanshift_3D(dx, dy, dz, xstart, ystart, zstart, delta_broad, maxiter);
rho = length(x1)^2 / area;
[~, ~, ~, ntruepairs_broad, ~, ~, ~, loc_error_auto, loc_error_auto_z] = meanshift_3D_error(dx, dy, dz, xshift_broadsweep(end), yshift_broadsweep(end), zshift_broadsweep(end), delta_broad, rmax, rho, numel(x1), numel(x2));
loc_prec = loc_error_auto/sqrt(2);
loc_prec_z = loc_error_auto_z/sqrt(2);
delta_narrow = loc_prec*specs.delta_narrow_ratio;
mean_nframeson = (ntruepairs_broad/2)/length(x1) + 1; % relevant to the error calculation, too

parfor i = 1:nTimeBin - 1
    refdata  = getnthdata(i);

    xshift_temp = zeros(nTimeBin, 1); 
    yshift_temp = zeros(nTimeBin, 1);
    zshift_temp = zeros(nTimeBin, 1);
    dxshift_temp = zeros(nTimeBin, 1); 
    dyshift_temp = zeros(nTimeBin, 1);
    dzshift_temp = zeros(nTimeBin, 1);
    
    ntruepairs_temp = zeros(nTimeBin, 1); nfalsepairs_temp = zeros(nTimeBin, 1);
    ntruepairs_z_temp = zeros(nTimeBin, 1); nfalsepairs_z_temp = zeros(nTimeBin, 1);
    loc_error_temp = zeros(nTimeBin, 1); loc_error_z_temp = zeros(nTimeBin, 1);
    iter_toconverge_temp = zeros(nTimeBin, 1);
    for j = i+1:nTimeBin
        thisdata = getnthdata(j);
        
        x1 = refdata(:,1); y1 = refdata(:,2); z1 = refdata(:,3); t1 = refdata(:,4);
        x2 = thisdata(:,1); y2 = thisdata(:,2); z2 = thisdata(:,3); t2 = thisdata(:,4);
        
        % translate this data by the previous shifts
        xoffset = xshift_temp(j-1); yoffset = yshift_temp(j-1); zoffset = zshift_temp(j-1);      
        x2 = x2 - xoffset; y2 = y2 - yoffset; z2 = z2 - zoffset;
        taumin = min(min(t1) - max(t2), min(t1) - max(t2)); taumax = -taumin;
        
        [dx, dy, dz] = crosspairs_3D(x1, y1, z1, t1, x2, y2, z2, t2, rmax, taumin, taumax, noutmax);
        
        xstart = 0; ystart = 0; zstart = 0;
        if broadsweep % necessary if you expect large shifts between consecutive bins
            [xshift_broadsweep, yshift_broadsweep, zshift_broadsweep] = meanshift_3D(dx, dy, dz, xstart, ystart, zstart, delta_broad, maxiter);
            xstart = xshift_broadsweep(end); ystart = yshift_broadsweep(end); zstart = zshift_broadsweep(end);
            if isnan(xstart)
                xstart = 0; ystart = 0; zstart = 0;
            end
        end
        
        [xshift_narrow, yshift_narrow, zshift_narrow] = meanshift_3D(dx, dy, dz, xstart, ystart, zstart, delta_narrow, maxiter);
        iter_toconverge_temp(j) = length(xshift_narrow) - 1;
        xshift_temp(j) = xshift_narrow(end) + xoffset;
        yshift_temp(j) = yshift_narrow(end) + yoffset;
        zshift_temp(j) = zshift_narrow(end) + zoffset;
        
        if isnan(xshift_temp(j))
            xshift_temp(j) = xoffset;
            yshift_temp(j) = yoffset;
            zshift_temp(j) = zoffset;
        end
        
        if calc_error
            rho = length(x1) * length(x2) / area;
            [dxshift_temp(j), dyshift_temp(j), dzshift_temp(j), ntruepairs_temp(j), nfalsepairs_temp(j), ntruepairs_z_temp(j), nfalsepairs_z_temp(j), loc_error_temp(j), loc_error_z_temp(j)] = meanshift_3D_error(dx, dy, dz, ...
                xshift_narrow(end), yshift_narrow(end), zshift_narrow(end), delta_narrow, rmax, rho, numel(x1), numel(x2), mean_nframeson, loc_prec, loc_prec_z);
            if isnan(dxshift_temp(j))
                dxshift_temp(j) = rmax;
                dyshift_temp(j) = rmax;
                dzshift_temp(j) = rmax;
            end
        end
    end
    
    % updating the matrices with the temporary rows
    xshift(i, :) = xshift_temp; yshift(i, :) = yshift_temp; zshift(i, :) = zshift_temp;
    dxshift(i, :) = dxshift_temp; dyshift(i, :) = dyshift_temp; dzshift(i, :) = dzshift_temp;
    ntruepairs(i, :) = ntruepairs_temp; nfalsepairs(i, :) = nfalsepairs_temp;
    ntruepairs_z(i, :) = ntruepairs_z_temp; nfalsepairs_z(i, :) = nfalsepairs_z_temp;
    loc_error(i, :) = loc_error_temp; loc_error_z(i, :) = loc_error_z_temp;
    iter_toconverge(i, :) = iter_toconverge_temp;   
%     fprintf('Finished with timebin %i\n', i)
end

% check if median predicted error is > loc_prec/4. If so, we may have chosen too
% few frames for alignment.

med_error = median(dxshift(dxshift > 0));
if med_error > loc_prec/4
    warning('Errors between alignments are large: consider decreasing the number of time bins.')
end

% redundant calculation is adapted from code in "Wang et al."

% total number of iterations is nTimeBin*(nTimeBin-1)/2
nelements = (nTimeBin^2-nTimeBin)/2;
A = zeros(nelements, 2);
R = zeros(nelements, 3);
w = zeros(nelements, 1);
count = 1;
for i=1:nTimeBin-1
    for j=i+1:nTimeBin
        A(count, i) = -1;
        A(count, j) = 1;
        R(count, 1) = xshift(i,j);
        R(count, 2) = yshift(i,j);
        R(count, 3) = zshift(i,j);
        if calc_error
            w(count, 1) = 1/dxshift(i,j)^2;
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
    rowerr(i,1) = sqrt(err(i,1)^2+err(i,2)^2+err(i,3)^2);
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

drift = vertcat([0,0,0], D);
stderr = vertcat([0,0,0], stderr);
xinterp = interp1(midtiming, drift(:,1), timings, interp_method, 'extrap');
yinterp = interp1(midtiming, drift(:,2), timings, interp_method, 'extrap');
zinterp = interp1(midtiming, drift(:,3), timings, interp_method, 'extrap');

% info and diagnostics
drift_info.xshift = xshift;
drift_info.yshift = yshift;
drift_info.zshift = zshift;
drift_info.dxshift = dxshift;
drift_info.dyshift = dyshift;
drift_info.dzshift = dzshift;
drift_info.xinterp = xinterp;
drift_info.yinterp = yinterp;
drift_info.zinterp = zinterp;
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
drift_info.ntruepairs_z = ntruepairs_z;
drift_info.nfalsepairs_z = nfalsepairs_z;
drift_info.loc_error = loc_error;
drift_info.loc_error_z = loc_error_z;
drift_info.iter_toconverge = iter_toconverge;

drift_info.delta_narrow = delta_narrow;
drift_info.loc_prec = loc_prec; 
drift_info.loc_prec_z = loc_prec_z; 
drift_info.med_error = med_error;

% Shift the data by the interpolated drift values at each time point.
shifteddata = olddata;
for i = 1:nframes
    index = find(data(:, 4) == i);
    shifteddata(index, 1) = olddata(index, 1) - xinterp(i);
    shifteddata(index, 2) = olddata(index, 2) - yinterp(i);
    shifteddata(index, 3) = olddata(index, 3) - zinterp(i);
end
shifteddata(:,4) = olddata(:,4);


