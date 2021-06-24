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

function [dxshift, dyshift, ntruepairs, nfalsepairs, loc_error] = meanshift_2D_error(dx, dy, xshift, yshift, delta, rmax, rho, mean_nframeson, loc_prec)
% Provides an error estimate for the 2D mean shift approach
dx_corr = dx - xshift;
dy_corr = dy - yshift;
r = sqrt(dx_corr.^2+dy_corr.^2);

ncandidatepairs = length(find(r < delta));

% choose a bin size in r for the correlation function.
if ncandidatepairs < 100
    rbinsize = delta/10;
else
    rbinsize = delta/20;
end

rbins = 0:rbinsize:rmax;
Nr = histcounts(r, rbins);
rmid = 0.5 * (rbins(1:end-1) + rbins(2:end));

inner_rads = rbins(1:end-1);
outer_rads = rbins(2:end); 
areas = pi*(outer_rads.^2 - inner_rads.^2);

Nr_exp = rho.*areas;
Nr_norm = Nr./Nr_exp;

deltabins = floor(delta / rbinsize);
baseline = median(Nr_norm);
ntruepairs = sum(Nr(1:deltabins)-baseline*Nr_exp(1:deltabins));
nfalsepairs = sum(Nr(1:deltabins)) - ntruepairs;
if ntruepairs < 0
    ntruepairs = 1;
end
if nfalsepairs < 0
    nfalsepairs = 0; 
end

% fit the correlation function to a Gaussian
fitgauss = fittype(@(A, s, c, r) A*exp(-r.^2/(2*s.^2)) + c,...
    'coefficients', {'A', 's', 'c'},...
    'indep', 'r', 'dep', 'z');
fgo = fitoptions(fitgauss);
fgo.StartPoint = [max(Nr_norm(1:deltabins)), delta/2, min(Nr_norm(1:deltabins))];
fgo.Lower = [0, 0, 0];
fgo.Upper = [Inf, rmax, Inf];

F = fit(rmid(1:deltabins)', Nr_norm(1:deltabins)', fitgauss, fgo);
% width of the correlation function gives the error for each pair of points
loc_error = F.s;
CI = confint(F, .68);
d = .5*(diff(CI, 1)); 

if isnan(d(2)) || d(2) > delta/2
    loc_error = delta;
end

if exist('loc_prec', 'var')
    loc_error = max(loc_error, sqrt(2)*loc_prec);
end

if numel(dx) == 1
    loc_error = NaN;
end

totalnpairs = ntruepairs + nfalsepairs;
if ~exist('mean_nframeson', 'var')
    mean_nframeson = 1;
end

dxshift = 1/sqrt(2)*sqrt(2*loc_error^2*ntruepairs + delta^2/2*nfalsepairs*mean_nframeson^2)/totalnpairs*(1+nfalsepairs/totalnpairs);
dyshift = dxshift;
