function nerve = buildQuicklookNerve(p, cuffSites) % Build a sampled nerve surface with local rigid-plane contact.
%% 1. Nerve samples and contact stations
R = p.geometry.radius_m; % Read the undeformed circular nerve radius.
L = p.domain.length_m; % Read the modelled length along the nerve.
depth = p.electrode.cuffIndentation_m; % Read the prescribed maximum flattening.
transition = p.electrode.cuffTransition_m; % Read the axial transition length outside the pad.
theta = linspace(0, 2*pi, 361); % Sample the external nerve surface around each section.
z_m = linspace(-L/2, L/2, 301); % Create regular longitudinal samples across the nerve.
normals = [cosd(cuffSites.theta_deg), sind(cuffSites.theta_deg)]; % Get outward directions from the device data.
assert(depth > 0 && depth < R && transition > 0, 'Invalid cuff flattening dimensions.'); % Reject invalid planes or tapers.
assert(all(abs(cuffSites.z_m) + cuffSites.sideZ_m/2 + transition <= L/2), ...
    'Cuff region lies beyond the nerve.'); % Keep all cuffs and their tapers within the modelled segment.
for k = 1:height(cuffSites) % Refine the axial surface where each pad and taper sit.
    span = cuffSites.sideZ_m(k)/2 + transition; % Include the pad width and its transition on both sides.
    zHere = cuffSites.z_m(k) + linspace(-span, span, 41); % Sample the entire local deformation.
    zEdges = cuffSites.z_m(k) + [-1, 0, 1]*cuffSites.sideZ_m(k)/2; % Include exact axial pad edges.
    z_m = [z_m, zHere, zEdges]; %#ok<AGROW> % Merge the local samples into the global grid.
end % Finish refining stations.
z_m = unique(z_m(:)); % Sort axial positions and remove duplicates.
directions = [cos(theta(:)), sin(theta(:))]; % Represent rays in the transverse plane.
radius_m = repmat(nerveBoundaryRadius(p,directions)',numel(z_m),1); % Start with the configured oval at every level.
%% 2. Flatten the nerve and place the rigid active surfaces
cuffFaces_m = zeros(4, 3, height(cuffSites)); % Allocate four vertices for each finite-area patch.
rPeri_m = max(p.geometry.fascicles_m(:,3:4),[],2) + p.geometry.tPeri_m; % Include the perineurium in containment checks.
centres_m = p.geometry.fascicles_m(:,1:2); % Keep the provisional fascicles straight and fixed.
for k = 1:height(cuffSites) % Clip the circular surface against each fixed cuff plane.
    n = normals(k,:); % Get the outward unit normal of this cuff.
    tangent = [-n(2), n(1)]; % Get the in-plane transverse direction of this cuff.
    t = max(0, min(1, (abs(z_m - cuffSites.z_m(k)) - cuffSites.sideZ_m(k)/2) / transition)); % Scale distance through the taper.
    support_m = nerveBoundaryRadius(p,n); % Find the undeformed outer radius in this pad direction.
    plane_m = support_m - depth*(1 - t.^2.*(3 - 2*t)); % Use a flat plateau and a smooth axial transition.
    projection = directions*n'; % Project each transverse ray onto the cuff normal.
    forward = projection > 0; % Only outward-facing rays can meet this plane.
    radius_m(:,forward) = min(radius_m(:,forward), plane_m ./ projection(forward)'); % Trim rays at the plane.
    halfWidth = cuffSites.sideXY_m(k)/2; % Read half the transverse pad width.
    halfLength = cuffSites.sideZ_m(k)/2; % Read half the axial pad length.
    assert(all(vecnorm((support_m-depth)*n + [-1;1]*halfWidth*tangent,2,2) <= max(p.geometry.nerveSemiAxes_m)), 'Cuff pad exceeds nerve bounds.'); % Reject grossly oversized pads.
    assert(all(centres_m*n' + rPeri_m < support_m-depth), 'Cuff plane intersects a fascicle.'); % Skip infeasible deformation inputs.
    xy = (support_m-depth)*n + [-1; 1; 1; -1]*halfWidth*tangent; % Put all four pad vertices on the flat surface.
    zz = cuffSites.z_m(k) + [-1; -1; 1; 1]*halfLength; % Set the pad's finite axial interval.
    cuffFaces_m(:,:,k) = [xy, zz]; % Store the active face without a bulk electrode.
end % Finish applying the cuff stations.
nerve.x_m = radius_m .* cos(theta); % Convert transverse radial samples to x.
nerve.y_m = radius_m .* sin(theta); % Convert transverse radial samples to y.
nerve.z_m = repmat(z_m, 1, numel(theta)); % Give every point its axial coordinate.
nerve.cuffFaces_m = cuffFaces_m; % Pass the rigid contacts to the renderer.
end % Return the nerve surface and finite-area cuff contacts.
