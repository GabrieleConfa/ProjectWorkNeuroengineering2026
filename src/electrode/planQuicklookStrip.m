function strip = planQuicklookStrip(p, sites, cuffSites, mode) % Place one transverse strip and its contact faces.
if mode=="deflect" && isfield(p.electrode,'routingModel') && p.electrode.routingModel=="beam" % Use the current physical model for new cases.
    strip=equilibriumStrip(p,sites,cuffSites); % Solve for the minimum-energy admissible shape.
    return % Older saved cases retain their original pad coordinates in the export.
end % Finish selecting the recorded model version.
%% 1. Device coordinates and collision envelope
[~, order] = sort(sites.r_m, 'descend'); % Identify an outer tutor contact to orient the strip axis.
u = sites.centre_m(order(1),1:2) / sites.r_m(order(1)); % Get the nominal transverse insertion direction.
v = [-u(2), u(1)]; % Define the lateral direction used for candidate deflections.
z = sites.z_m(1); % Keep the station at its measured longitudinal coordinate.
R = max(p.geometry.nerveSemiAxes_m); % Use the outer scale for the deflection family.
entryRadius_m = nerveBoundaryRadius(p,-u); % Find the configured nerve edge at the entry side.
exitRadius_m = nerveBoundaryRadius(p,u); % Find the configured nerve edge at the exit side.
s_m = linspace(-entryRadius_m, exitRadius_m, 301)'; % Sample between the true oval surface intersections.
span_m = (entryRadius_m+exitRadius_m)/2; % Set the half-span for a dimensionless bending profile.
midpoint_m = (exitRadius_m-entryRadius_m)/2; % Locate the midpoint between the two boundary points.
siteS_m = sites.centre_m(:,1:2)*u'; % Preserve each tutor contact's nominal along-strip coordinate.
if ~isfield(p.electrode,'tipOverhang_m'), p.electrode.tipOverhang_m=0.05e-3; end % Support previously saved personal cases.
tip = max(siteS_m + sites.sideXY_m/2) + p.electrode.tipOverhang_m; % Extend beyond the furthest nominal pad edge.
assert(p.electrode.tipOverhang_m>0 && tip<exitRadius_m,'Lead tip must end inside the nerve beyond the last pad.'); % Reject incompatible device dimensions.
s_m = [s_m(s_m<tip);tip]; % Sample only the single continuous insertion from the entry wall to the interior tip.
fasc = p.geometry.fascicles_m; % Read the fixed provisional fascicles.
peri_m = max(fasc(:,3:4),[],2) + p.geometry.tPeri_m; % Use the outer perineurium as the fascicle boundary.
margin_m = p.electrode.stripWidth_m/2; % Derive the geometric collision envelope from the current strip width.
%% 2. Search for a smooth lateral route or retain a straight crossing
offsets_m = [0, p.electrode.maxDeflection_m * reshape([1:160; -(1:160)], 1, []) / 160]; % Try small shifts first, in either direction.
skew=zeros(size(offsets_m)); wave=skew; % Retain the original candidate family for archived solutions.
if mode=="deflect" && isfield(p.electrode,'routeSeed') % Sample independent smooth insertion scenarios for new cases.
    stream=RandStream('mt19937ar','Seed',p.electrode.routeSeed+100+sites.section(1)); % Reproduce this axial station irrespective of call order.
    count=800; % Bound the number of proposed geometries before the deterministic fallback.
    amplitude=p.electrode.maxDeflection_m*(2*rand(stream,1,count)-1); % Sample both sides of a possible obstacle.
    offsets_m=[0,amplitude,offsets_m(2:end)]; % Retain a straight route whenever the anatomy permits it.
    skew=[0,1.5*(2*rand(stream,1,count)-1),skew(2:end)]; % Vary the location of peak bending.
    wave=[0,.65*(2*rand(stream,1,count)-1),wave(2:end)]; % Allow asymmetric bends at successive encounters.
end % Every sampled proposal must pass the same tissue and device checks.
if mode == "penetrate" % Select the intentional straight-line scenario.
    offsets_m = 0; % Keep the tutor's nominal path and report intersections.
end % Finish selecting the route search.
for attempt=1:numel(offsets_m) % Test the candidate trajectories in their reproducible order.
    offset_m=offsets_m(attempt); % Read the signed displacement for this trajectory.
    q = (s_m-midpoint_m)/span_m; % Normalize the insertion position between the two surfaces.
    shape = 1+skew(attempt)*q+wave(attempt)*sin(2*pi*q); % Vary local curvature across the insertion.
    bump = (1 - q.^2).^2.*shape; % Keep the strip entrances and their directions fixed.
    path_m = s_m*u + offset_m*bump*v; % Build the candidate centreline in the transverse plane.
    if max(abs(offset_m*bump))>p.electrode.maxDeflection_m+eps, continue; end % Enforce the editable lateral travel limit.
    distance_m = segmentDistances(path_m(1:end-1,:), path_m(2:end,:), fasc(:,1:2)); % Check complete polyline segments.
    padXY = zeros(2, 2, height(sites)); % Allocate the two transverse endpoints of each planar pad.
    padDistance_m = zeros(height(sites), size(fasc,1)); % Allocate contact-to-fascicle distances.
    for k = 1:height(sites) % Construct each finite pad on the deformed path.
        s = siteS_m(k); % Read this contact's fixed station along the nominal insertion line.
        xy = interp1(s_m, path_m, s, 'linear'); % Locate the pad centre on the sampled strip.
        qPad = (s-midpoint_m)/span_m; % Normalize this pad position along the insertion.
        slope = offset_m/span_m*((-4*qPad*(1-qPad^2))*(1+skew(attempt)*qPad+wave(attempt)*sin(2*pi*qPad))+(1-qPad^2)^2*(skew(attempt)+2*pi*wave(attempt)*cos(2*pi*qPad))); % Differentiate the lateral displacement at the pad.
        tangent = u + slope*v; % Match the pad direction to the strip tangent.
        tangent = tangent / norm(tangent); % Preserve the original transverse pad width.
        padXY(:,:,k) = xy + [-1; 1]*sites.sideXY_m(k)/2*tangent; % Build the straight edges of a flat contact.
        padDistance_m(k,:) = segmentDistances(padXY(1,:,k), padXY(2,:,k), fasc(:,1:2)); % Test the full contact width.
    end % Finish checking the pads.
    pathRadius_m = cuffDeformedRadius(p,cuffSites,z,path_m ./ max(vecnorm(path_m,2,2),eps)); % Check the actual locally indented contour.
    inside = all(vecnorm(path_m,2,2) <= pathRadius_m' + 1e-12); % Keep the strip inside the nerve at its own section.
    padPoints = reshape(permute(padXY,[1 3 2]), [], 2); % Collect all transverse pad endpoints.
    padRadius_m = cuffDeformedRadius(p,cuffSites,z,padPoints ./ max(vecnorm(padPoints,2,2),eps)); % Check the deformed contour at each finite pad edge.
    padInside = all(vecnorm(padPoints,2,2) <= padRadius_m' + 1e-12); % Keep all active faces within the nerve.
    hasClearance = all(distance_m > peri_m' + margin_m, 'all') && all(padDistance_m > peri_m' + margin_m, 'all'); % Require clearance from every fascicle.
    if inside && padInside && (hasClearance || mode == "penetrate") % Accept only an admissible scenario.
        break % Keep the first route passing the checks.
    end % Move to the next offset when this candidate fails.
end % Finish searching for a route.
assert(inside && padInside && (hasClearance || mode == "penetrate"), ...
    'Section %+d has no route within the configured deflection family.', sites.section(1)); % Report an infeasible station explicitly.
%% 3. Export the path, pad faces and intersected fascicles
strip.path_m = [path_m, repmat(z, size(path_m,1), 1)]; % Store the 3D centreline.
strip.axisXY = u; % Keep the nominal entry direction for straight-reference plots.
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
