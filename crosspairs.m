function [dxout, dyout, dtout, status] = crosspairs(x1, y1, t1, x2, y2, t2, rmax,...
        taumin, taumax, noutmax);

% Sort based on x
[x1, s1] = sort(x1(:));
x1 = x1';
[x2, s2] = sort(x2(:));
x2 = x2';

y1 = y1(s1);
t1 = t1(s1);
y2 = y2(s2);
t2 = t2(s2);

[dxout, dyout, dtout, status] = Fcrosspairs(x1, y1, t1, x2, y2, t2,...
                                    rmax, taumin, taumax, int64(noutmax));