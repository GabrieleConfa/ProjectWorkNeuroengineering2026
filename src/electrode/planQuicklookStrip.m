function strip = planQuicklookStrip(p, sites, cuffSites, mode) % Place one transverse strip and its contact faces.
%% 1. Device coordinates and collision envelope
[~, order] = sort(sites.r_m, 'descend'); % Identify an outer tutor contact to orient the strip axis.
u = sites.centre_m(order(1),1:2) / sites.r_m(order(1)); % Get the nominal transverse insertion direction.
v = [-u(2), u(1)]; % Define the lateral direction used for candidate deflections.
z = sites.z_m(1); % Keep the station at its measured longitudinal coordinate.
R = p.geometry.radius_m; % Read the circular nerve radius away from cuff clipping.
s_m = linspace(-R, R, 301)'; % Sample the strip from one nerve surface to the opposite surface.
siteS_m = sites.centre_m(:,1:2)*u'; % Preserve each tutor contact's nominal along-strip coordinate.
fasc = p.geometry.fascicles_m; % Read the fixed provisional fascicles.
peri_m = fasc(:,3) + p.geometry.tPeri_m; % Use the outer perineurium as the fascicle boundary.
margin_m = p.electrode.stripHalfEnvelope_m; % Inflate obstacles by the illustrative strip half-envelope.
cuffN = [cosd(cuffSites.theta_deg), sind(cuffSites.theta_deg)]; % Read the station's outward cuff normals.
cuffPlane_m = R - p.electrode.cuffIndentation_m; % Use the same cuff plane as the nerve surface.
%% 2. Search for a smooth lateral route or retain a straight crossing
offsets_m = [0, p.electrode.maxDeflection_m * reshape([1:160; -(1:160)], 1, []) / 160]; % Try small shifts first, in either direction.
if mode == "penetrate" % Select the intentional straight-line scenario.
    offsets_m = 0; % Keep the tutor's nominal path and report intersections.
end % Finish selecting the route search.
for offset_m = offsets_m % Evaluate offsets in order of magnitude.
    bump = (1 - (s_m/R).^2).^2; % Keep the strip entrances and their directions fixed.
    path_m = s_m*u + offset_m*bump*v; % Build the candidate centreline in the transverse plane.
    distance_m = segmentDistances(path_m(1:end-1,:), path_m(2:end,:), fasc(:,1:2)); % Check complete polyline segments.
    padXY = zeros(2, 2, height(sites)); % Allocate the two transverse endpoints of each planar pad.
    padDistance_m = zeros(height(sites), size(fasc,1)); % Allocate contact-to-fascicle distances.
    for k = 1:height(sites) % Construct each finite pad on the deformed path.
        s = siteS_m(k); % Read this contact's fixed station along the nominal insertion line.
        xy = interp1(s_m, path_m, s, 'linear'); % Locate the pad centre on the sampled strip.
        slope = -4*offset_m*s/R^2 * (1 - (s/R)^2); % Differentiate the lateral displacement at the pad.
        tangent = u + slope*v; % Match the pad direction to the strip tangent.
        tangent = tangent / norm(tangent); % Preserve the original transverse pad width.
        padXY(:,:,k) = xy + [-1; 1]*sites.sideXY_m(k)/2*tangent; % Build the straight edges of a flat contact.
        padDistance_m(k,:) = segmentDistances(padXY(1,:,k), padXY(2,:,k), fasc(:,1:2)); % Test the full contact width.
    end % Finish checking the pads.
    inside = all(vecnorm(path_m,2,2) <= R + 1e-12) && all(path_m*cuffN' <= cuffPlane_m, 'all'); % Keep the intraneural centreline in the flattened station.
    padPoints = reshape(permute(padXY,[1 3 2]), [], 2); % Collect all transverse pad endpoints.
    padInside = all(vecnorm(padPoints,2,2) <= R) && all(padPoints*cuffN' <= cuffPlane_m, 'all'); % Keep finite contacts inside the tissue.
    hasClearance = all(distance_m > peri_m' + margin_m, 'all') && all(padDistance_m > peri_m' + margin_m, 'all'); % Require clearance from every fascicle.
    if inside && padInside && (hasClearance || mode == "penetrate") % Accept only an admissible scenario.
        break % Keep the first route passing the checks.
    end % Move to the next offset when this candidate fails.
end % Finish searching for a route.
assert(inside && padInside && (hasClearance || mode == "penetrate"), ...
    'Section %+d has no route within the configured deflection family.', sites.section(1)); % Report an infeasible station explicitly.
%% 3. Export the path, pad faces and intersected fascicles
strip.path_m = [path_m, repmat(z, size(path_m,1), 1)]; % Store the 3D centreline.
strip.section = sites.section(1); % Record the tutor's section identifier.
strip.offset_m = offset_m; % Record the chosen lateral deflection.
strip.envelopeFascicles = find(any([distance_m; padDistance_m] <= peri_m' + margin_m, 1)); % List geometric overlaps with the provisional strip envelope.
strip.padFaces_m = zeros(4, 3, height(sites)); % Allocate the 3D pad vertices.
for k = 1:height(sites) % Extend each transverse contact over its axial length.
    zz = z + [-1; -1; 1; 1]*sites.sideZ_m(k)/2; % Read pad edges along the nerve axis.
    strip.padFaces_m(:,:,k) = [padXY([1,2,2,1],:,k), zz]; % Export a single ideal contact surface per pad.
end % Finish constructing pad surfaces.
end % Return one completed station.
function distance_m = segmentDistances(a, b, centres) % Measure the shortest distances from segments to points.
d = b - a; % Set the vector along each segment.
offsetX = centres(:,1)' - a(:,1); % Set the x component from each segment origin to each centre.
offsetY = centres(:,2)' - a(:,2); % Set the y component from each segment origin to each centre.
t = (offsetX.*d(:,1) + offsetY.*d(:,2)) ./ sum(d.^2,2); % Project centres onto the infinite lines.
t = max(0, min(1, t)); % Clamp the closest point to the segment endpoints.
distance_m = hypot(offsetX-t.*d(:,1), offsetY-t.*d(:,2)); % Measure exact planar clearance for the polyline and pads.
end % Return one distance per segment and fascicle.
