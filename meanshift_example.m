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

%% Load in the data
close all
data = csvread('Data/NPCdata_meanshift_example.csv');
index1 = find(1 <= data(:, 3) &  data(:, 3) <= 500);
data1 = data(index1, :);
index2 = find(10001 <= data(:, 3) &  data(:, 3) <= 10500);
data2 = data(index2, :);

x1 = data1(:,1); y1 = data1(:,2); t1 = data1(:,3);
x2 = data2(:,1); y2 = data2(:,2); t2 = data2(:,3);

%% Extract pairwise displacements between the two datasets
taumin = 0; taumax = max(t2)-min(t1); rmax = 500; noutmax = numel(x1)*numel(x2);
[dx, dy] = crosspairs(x1, y1, t1, x2, y2, t2, rmax, taumin, taumax, noutmax);

%% Computing the Shift Between the datasets
% Initially choose a large delta due to the large displacement in y.
delta_broad = 150; 
xstart = 0; ystart = 0; maxiter = 200;
[xshift_broad_full, yshift_broad_full] = meanshift_2D(dx, dy, xstart, ystart, delta_broad, maxiter);
xshift_broad = xshift_broad_full(end); yshift_broad = yshift_broad_full(end); 

%% Choosing More Precise Mean Shift Parameters
% Now that we have identified the peaks, we can choose a narrower delta for
% a more precise estimate.
delta_narrow = 50;
[xshift_full, yshift_full] = meanshift_2D(dx, dy, xstart, ystart, delta_narrow, maxiter);
xshift = xshift_full(end); yshift = yshift_full(end); 

% Visualize the displacements and the computed shift.
figure; h1 = scatter(dx, dy, 1, 'filled', 'MarkerFaceAlpha', 0.1, 'MarkerEdgeAlpha', 0.1);
hold on; h2 = plot(0, 0, 'k+', 'MarkerSize', 20);
hold on; h3 = plot(xshift, yshift, 'k*', 'MarkerSize', 20);
xlabel('dx (nm)')
ylabel('dy (nm)')
axis equal
lgd = legend([h2 h3],{'Starting Point','Calculated Shift'});
lgd.Position = [0.6479    0.7571    0.2357    0.1320];

% Compute error estimates
rho = numel(x1)*numel(x2) / ((max([x1; x2])-min([x1; x2]))*(max([y1; y2])-min([y1; y2])));
[dxshift, dyshift] = meanshift_2D_error(dx, dy, xshift, yshift, delta_narrow, rmax, rho);
