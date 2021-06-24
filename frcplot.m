% Copyright (C) 2021 Frank Fazekas
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

% (C) Copyright 2012               Quantitative Imaging Group
%     All rights reserved          Faculty of Applied Physics
%                                  Delft University of Technology
%                                  Lorentzweg 1
%                                  2628 CJ Delft
%                                  The Netherlands
% Robert Nieuwenhuizen, Oct 2012

function frcplot(fire_value, fireL, fireH, pixelsize, superzoom, frc_curve)

fprintf('FIRE value %2.1f +- %2.2f [px]\n', fire_value, (fireL-fireH)/2);
fprintf('FIRE value %2.1f +- %2.2f [nm]\n', fire_value*pixelsize/superzoom, (fireL-fireH)/2*pixelsize/superzoom);
figure;
qmax = 0.5/(pixelsize/superzoom);
plot(linspace(0,qmax*sqrt(2), length(frc_curve)), frc_curve,'-')
q = linspace(0,qmax*sqrt(2), length(frc_curve));
xlim([0,qmax])
hold on
plot([0 qmax],[1/7 1/7],'r-');
plot([0 qmax],[0 0],'k--'); hold off
xlabel('spatial frequency (nm^{-1})')
ylabel('FRC')
title('Fourier Ring Correlation curve')