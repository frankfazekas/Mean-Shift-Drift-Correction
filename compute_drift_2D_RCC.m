function [shifteddata, drift_info] = compute_drift_2D_RCC(data, specs, timings)
% adapted from supplementary software in "Wang et al."

addpath('RCC')
 
% sort by frame
olddata = data;
data = sort(data, 3);
nframes = max(data(:, 3));

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

% initialize variables
xshift = zeros(nTimeBin, nTimeBin);
yshift = zeros(nTimeBin, nTimeBin);

imsize = max(max(data(:,1)), max(data(:,2)));

% can use for loop instead of parfor if no parallel capabilities are available
parfor i = 1:nTimeBin - 1  
    index = data(:,3)>=firstframe(i) & data(:,3)<=lastframe(i);
    imorigin = BinLocalizations(data(index,1:2), imsize, zoomfactor);
    
    autocorr = CrossCorrelation(imorigin,imorigin);
    [yorigin, xorigin] = GaussianFit(autocorr);
    
    xshift_temp = zeros(nTimeBin, 1); 
    yshift_temp = zeros(nTimeBin, 1);
    for j = i+1:nTimeBin
        index = data(:,3)>=firstframe(j) & data(:,3)<=lastframe(j);
        imbin = BinLocalizations(data(index,1:2), imsize, zoomfactor);
        
        corr = CrossCorrelation(imorigin,imbin);
        [ycross,xcross] = GaussianFit(corr);
        
        xshift_temp(j) = (xorigin - xcross)/zoomfactor;
        yshift_temp(j) = (yorigin - ycross)/zoomfactor; 
          
    end
    
    % updating the matrices with the temporary rows
    xshift(i, :) = xshift_temp; yshift(i, :) = yshift_temp;
%     fprintf('Finished with time bin %i\n', i)
end

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

% Shift the data by the interpolated drift values at each time point.
shifteddata = olddata;
for i = 1:nframes
    index = find(data(:, 3) == i);
    shifteddata(index, 1) = olddata(index, 1) - xinterp(i);
    shifteddata(index, 2) = olddata(index, 2) - yinterp(i);
end
shifteddata(:,3) = olddata(:,3);