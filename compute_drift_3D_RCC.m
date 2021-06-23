function [shifteddata, drift_info] = compute_drift_3D_RCC(data, specs, timings)
% adapted from supplementary software in "Wang et al."

addpath('RCC')

% sort by frame

olddata = data;
data = sort(data, 4);
nframes = max(data(:, 4));

% Get parameters
nTimeBin = round(specs.nTimeBin);
binwidth = floor(nframes/nTimeBin);
psize = specs.psize;
zoomfactor = 1/psize;
outlier_error = specs.outlier_error;
interp_method = specs.interp_method;

% Compute which frames belong to which time bins
binspacing = (nframes - binwidth)/(nTimeBin-1);
firstframe = round(1 + (0:nTimeBin-1)*binspacing);
lastframe = min(nframes, firstframe + binwidth - 1);

xshift = zeros(nTimeBin, nTimeBin);
yshift = zeros(nTimeBin, nTimeBin);
zshift_xproj = zeros(nTimeBin, nTimeBin);
zshift_yproj = zeros(nTimeBin, nTimeBin);

imsize = max(max(data(:,1)), max(data(:,2)));
orig_zmin = min(data(:,3));
data(:,3) = data(:,3) - orig_zmin;
% can use for loop instead of parfor if no parallel capabilities are available
parfor i = 1:nTimeBin - 1
    index = data(:,4)>=firstframe(i) & data(:,4)<=lastframe(i);
    imorigin = BinLocalizations(data(index,1:2), imsize, zoomfactor);
    imorigin_xproj = BinLocalizations(data(index,2:3), imsize, zoomfactor);
    imorigin_yproj = BinLocalizations(data(index,[1 3]), imsize, zoomfactor);
    
    autocorr = CrossCorrelation(imorigin,imorigin);
    [yorigin, xorigin] = GaussianFit(autocorr);
    autocorr_xproj = CrossCorrelation(imorigin_xproj,imorigin_xproj);
    [zorigin_xproj, yorigin_xproj] = GaussianFit(autocorr_xproj);
    autocorr_yproj = CrossCorrelation(imorigin_yproj,imorigin_yproj);
    [zorigin_yproj, xorigin_yproj] = GaussianFit(autocorr_yproj);
    
    xshift_temp = zeros(nTimeBin, 1); 
    yshift_temp = zeros(nTimeBin, 1);
    zshift_temp_xproj = zeros(nTimeBin, 1);
    zshift_temp_yproj = zeros(nTimeBin, 1);
    for j = i+1:nTimeBin
        index = data(:,4)>=firstframe(j) & data(:,4)<=lastframe(j);
        
        imbin = BinLocalizations(data(index,1:2), imsize, zoomfactor);
        
        corr = CrossCorrelation(imorigin,imbin);
        [ycross,xcross] = GaussianFit(corr);
        
        xshift_temp(j) = (xorigin - xcross)/zoomfactor;
        yshift_temp(j) = (yorigin - ycross)/zoomfactor; 
        
        % x-projection to calculate z-drift
        imbin_xproj = BinLocalizations(data(index,2:3), imsize, zoomfactor);
        
        corr_xproj = CrossCorrelation(imorigin_xproj,imbin_xproj);
        [zcross_xproj,ycross_xproj] = GaussianFit(corr_xproj);
        
        zshift_temp_xproj(j) = (zorigin_xproj - zcross_xproj)/zoomfactor; 
        
        % y-projection to calculate z-drift
        imbin_yproj = BinLocalizations(data(index,[1 3]), imsize, zoomfactor);
        
        corr_yproj = CrossCorrelation(imorigin_yproj,imbin_yproj);
        [zcross_yproj,xcross_yproj] = GaussianFit(corr_yproj);
        
        zshift_temp_yproj(j) = (zorigin_yproj - zcross_yproj)/zoomfactor;     
        
    end
    % updating the matrices with the temporary rows
    xshift(i, :) = xshift_temp; yshift(i, :) = yshift_temp;
    zshift_xproj(i, :) = zshift_temp_xproj; zshift_yproj(i, :) = zshift_temp_yproj;   
%     fprintf('Finished with timebin %i\n', i)
end

zshift = (zshift_xproj + zshift_yproj)/2;

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
        w(count, 1) = 1;
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
drift_info.zshift_xproj = zshift_xproj;
drift_info.zshift_yproj = zshift_yproj;
drift_info.orig_zmin = orig_zmin;
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

% Shift the data by the interpolated drift values at each time point.
shifteddata = olddata;
for i = 1:nframes
    index = find(data(:, 4) == i);
    shifteddata(index, 1) = olddata(index, 1) - xinterp(i);
    shifteddata(index, 2) = olddata(index, 2) - yinterp(i);
    shifteddata(index, 3) = olddata(index, 3) - zinterp(i);
end
shifteddata(:,4) = olddata(:,4);


