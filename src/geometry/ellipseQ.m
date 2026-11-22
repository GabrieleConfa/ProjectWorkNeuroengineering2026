function Q = ellipseQ(xy_m, semiAxes_m, rotation_deg)
%ELLIPSEQ Quadratic form of a centred, rotated ellipse: Q < 1 inside, Q = 1 on it, Q > 1 outside.
%   Q = ellipseQ(xy_m, semiAxes_m, rotation_deg)
%
%   The nerve cross-section is an ellipse centred on the axis (D-10). In the frame of the ellipse
%   (x' along the first semi-axis) Q = (x'/a)^2 + (y'/b)^2. sqrt(Q) is the "elliptical radius": the
%   points with sqrt(Q) = s form the same ellipse scaled by s. A circle of radius R is a = b = R.
%
%   Inputs
%     xy_m          (n,2) points [x y] [m]
%     semiAxes_m    [a b] [m], a along the direction rotation_deg
%     rotation_deg  angle of the first semi-axis from x, counterclockwise [deg]
%   Output
%     Q             (n,1) [-]

c = cosd(rotation_deg);   s = sind(rotation_deg);
xl =  c * xy_m(:,1) + s * xy_m(:,2);                 % coordinates in the ellipse frame
yl = -s * xy_m(:,1) + c * xy_m(:,2);
Q = (xl / semiAxes_m(1)).^2 + (yl / semiAxes_m(2)).^2;
end
