function [plane_m, indentation_m] = cuffContactPlane(p, site) % Fit one rigid flat pad to the nerve's native outline.
%% 1. Resolve the rotated oval in the pad's local directions
n = [cosd(site.theta_deg),sind(site.theta_deg)]; % Read the outward normal from the supplied contact layout.
t = [-n(2),n(1)]; % Define the transverse direction across the active pad.
a = p.geometry.nerveSemiAxes_m(1); % Read the first anatomical semi-axis.
b = p.geometry.nerveSemiAxes_m(2); % Read the second anatomical semi-axis.
c = cosd(p.geometry.nerveRotation_deg); % Read the nerve's cross-section rotation.
s = sind(p.geometry.nerveRotation_deg); % Read the corresponding sine.
rotation = [c,s;-s,c]; % Transform global transverse points into oval coordinates.
nLocal = rotation*n'; % Express the pad normal in oval coordinates.
tLocal = rotation*t'; % Express the transverse pad direction in oval coordinates.
A = (nLocal(1)/a)^2+(nLocal(2)/b)^2; % Compute the quadratic coefficient along the normal.
B = nLocal(1)*tLocal(1)/a^2+nLocal(2)*tLocal(2)/b^2; % Compute the normal-lateral coupling.
C = (tLocal(1)/a)^2+(tLocal(2)/b)^2; % Compute the quadratic coefficient across the pad.
%% 2. Locate a rigid plane under both transverse pad edges
lateral_m = site.sideXY_m/2*[-1,1]; % Evaluate the two measured active-face edges.
discriminant = B^2*lateral_m.^2+A*(1-C*lateral_m.^2); % Solve the line-oval intersections analytically.
assert(all(discriminant>=0), 'Cuff pad is wider than the nerve at this location.'); % Reject an impossible contact footprint.
edgeNormal_m = (-B*lateral_m+sqrt(discriminant))/A; % Find the exterior nerve point at each lateral edge.
nativeCentre_m = 1/sqrt(A); % Find the original surface below the centre of the pad.
plane_m = min(edgeNormal_m)-p.electrode.cuffExtraIndentation_m; % Indent enough for the whole rigid pad to fit.
indentation_m = nativeCentre_m-plane_m; % Return the resulting central flattening depth.
assert(plane_m>0 && indentation_m>=0, 'Cuff plane must lie inside the original nerve contour.'); % Reject an invalid fitted plane.
end % Return the geometric contact plane and its indentation.
