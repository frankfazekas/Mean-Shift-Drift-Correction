function [dxout, dyout, dzout, dtout] = crosspairs_3D(x1, y1, z1, t1, x2, y2, z2, t2, rmax,...
        taumin, taumax, noutmax);

% Sort based on x1
[x1, s1] = sort(x1(:));
x1 = x1';
[x2, s2] = sort(x2(:));
x2 = x2';

y1 = y1(s1);
z1 = z1(s1);
t1 = t1(s1);
y2 = y2(s2);
z2 = z2(s2);
t2 = t2(s2);

[dxout, dyout, dzout, dtout] = Fcrosspairs_3D(x1, y1, z1, t1, x2, y2, z2, t2,...
                                    rmax, taumin, taumax, int64(noutmax));