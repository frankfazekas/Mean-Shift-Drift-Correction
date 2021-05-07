function [dxshift, dyshift, dzshift, ntruepairs, nfalsepairs, zntruepairs, znfalsepairs, loc_error, zloc_error] = meanshift_3D_error(dx, dy, dz, xshift, yshift, zshift, delta, rmax, rho, ndata1, ndata2, mean_nframeson, loc_prec, zloc_prec)
% Provides an error estimate for the 2D mean shift approach
closei_forr = find(abs(dz - zshift(end)) < delta);

dx_corr = dx(closei_forr) - xshift;
dy_corr = dy(closei_forr) - yshift;
r = sqrt(dx_corr.^2+dy_corr.^2);

% choose a bin size in r for the correlation function.
ncandidatepairs = length(find(r < delta));
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

if ~exist('mean_nframeson', 'var')
    mean_nframeson = 1;
end

falsepairsSD = delta/sqrt(2);

if exist('mean_nframeson', 'var')
    falsepairsSE = falsepairsSD/sqrt(nfalsepairs/(mean_nframeson)^2);       
else
    falsepairsSE = falsepairsSD/sqrt(nfalsepairs);
end 

if nfalsepairs <= 0
    falsepairsSE = 0;
end

xloc_trueerror = loc_error/sqrt(ntruepairs);
yloc_trueerror = loc_error/sqrt(ntruepairs);

totalnpairs = ntruepairs + nfalsepairs;
totalSE = sqrt((ntruepairs*xloc_trueerror)^2+(nfalsepairs*falsepairsSE)^2/2)/totalnpairs;
dxshift = totalSE * totalnpairs/ntruepairs;
dyshift = dxshift;

% Create correlation function in z.
closei_forz = find(r < delta);
dz_corr = dz(closei_forz) - zshift(end);
z = abs(dz_corr);
if numel(closei_forz) < 100
    zbinsize = 20;
elseif numel(closei_forz) < 1000
    zbinsize = 10;
else    
    zbinsize = 5;
end

zbins = 0:zbinsize:rmax;
[Nz, zedges] = histcounts(z, zbins);

zmid = 0.5 * (zedges(1:end-1) + zedges(2:end));
depth = max(dz) - min(dz);
zrho = ndata1 * ndata2 / depth;

if zrho == 0
    zntruepairs = 0;
    znfalsepairs = NaN;
    dzshift = NaN;
    zloc_error = delta;
    return
else
    Nz_exp = zrho*zbinsize/length(zbins);
end

try
    Nz_norm = Nz/Nz_exp;
catch
    Nz_norm = Nz;
end

zdeltabins = floor(2*delta / zbinsize);
zbaseline = median(Nz_norm(:));

try
    zntruepairs = sum(Nz(1:zdeltabins)-zbaseline*Nz_exp);
catch
    zntruepairs = 0;
    znfalsepairs = NaN;
    dzshift = NaN;
    zloc_error = delta;
    return
end

if zntruepairs <= 0
    zntruepairs = 0;
    znfalsepairs = NaN;
    dzshift = NaN;
    zloc_error = delta;
    return
end

% fit correlation function to Gaussian to get localization error
fitgaussz = fittype(@(A, s, c, z) A*exp(-z.^2/(2*s.^2)) + c,...
    'coefficients', {'A', 's', 'c'},...
    'indep', 'z', 'dep', 'a');

fgoz = fitoptions(fitgaussz);
fgoz.StartPoint = [max(Nz_norm(1:zdeltabins)), delta, median(Nz_norm(:))];
fgoz.Lower = [0, 0, 0];
fgoz.Upper = [Inf, 2*delta, Inf];

F_z = fit(zmid(1:zdeltabins)', Nz_norm(1:zdeltabins)', fitgaussz, fgoz);
zloc_error = F_z.s;
    
CI_z = confint(F_z, .68);
d_z = .5*(diff(CI_z, 1)); 

if isnan(d_z(2)) || d_z(2) > delta/2
    zloc_error = delta;
end

if exist('zloc_prec', 'var')
    zloc_error = max(zloc_error, sqrt(2)*zloc_prec);
end

if numel(dx) == 1
    zloc_error = NaN;
end

znfalsepairs = sum(Nz(1:zdeltabins)) - zntruepairs;
if znfalsepairs < 1
    znfalsepairs = 1;
end
if znfalsepairs < 0
    znfalsepairs = 0;
end

zfalsepairsSD = delta/sqrt(3);

if exist('auto_ratio', 'var')
    zfalsepairsSE = zfalsepairsSD/sqrt(znfalsepairs);        
else
    zfalsepairsSE = zfalsepairsSD/sqrt(znfalsepairs);
end 

if znfalsepairs == 0
    zfalsepairsSE = 0;
end

zloc_trueerror = zloc_error/sqrt(zntruepairs);
ztotalnpairs = zntruepairs + znfalsepairs;
ztotalSE = sqrt((zntruepairs*zloc_trueerror)^2+(znfalsepairs*zfalsepairsSE)^2)/ztotalnpairs;
dzshift = ztotalSE * ztotalnpairs/zntruepairs;
