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