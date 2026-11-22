function [xy_m, nOut] = nerveOutline(semiAxes_m, rotation_deg, nPoints)
%NERVEOUTLINE Points of the elliptical nerve boundary and the outward unit normals there.
%   [xy_m, nOut] = nerveOutline(semiAxes_m, rotation_deg, nPoints)
%
%   Used to draw the nerve and, densely sampled, to measure the distance of a point from the nerve
%   surface (fascicle generator, tests). Points are equally spaced in the ellipse parameter t,
%   x' = a cos t, y' = b sin t, then rotated; the outward normal is the gradient of the quadratic form,
%   (x'/a^2, y'/b^2), normalised.
%
%   Inputs
%     semiAxes_m    [a b] [m]; a = b gives a circle
%     rotation_deg  angle of the first semi-axis from x [deg]
%     nPoints       number of points (closed curve, the first point is not repeated)
%   Outputs
%     xy_m          (nPoints,2) boundary points [m]
%     nOut          (nPoints,2) outward unit normals

t = 2*pi * (0:nPoints-1)' / nPoints;
a = semiAxes_m(1);   b = semiAxes_m(2);
Rot = [cosd(rotation_deg) sind(rotation_deg); -sind(rotation_deg) cosd(rotation_deg)];   % row vectors
xy_m = [a * cos(t), b * sin(t)] * Rot;
nOut = [cos(t) / a, sin(t) / b] * Rot;
nOut = nOut ./ vecnorm(nOut, 2, 2);
end
