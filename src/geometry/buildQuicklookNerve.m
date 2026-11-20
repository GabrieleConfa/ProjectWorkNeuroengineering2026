function nerve = buildQuicklookNerve(p, cuffSites) % Build a sampled nerve surface with local rigid-plane contact.
%% 1. Nerve samples and contact stations
L = p.domain.length_m; % Read the modelled length along the nerve.
transition = p.electrode.cuffTransition_m; % Read the axial transition length outside the pad.
theta = linspace(0, 2*pi, 361); % Sample the external nerve surface around each section.
z_m = linspace(-L/2, L/2, 301); % Create regular longitudinal samples across the nerve.
normals = [cosd(cuffSites.theta_deg), sind(cuffSites.theta_deg)]; % Get outward directions from the device data.
assert(transition>0 && p.electrode.cuffLateralTransition_m>0, 'Cuff transitions must be positive.'); % Reject invalid shoulder widths.
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
radius_m = cuffDeformedRadius(p,cuffSites,z_m,directions); % Apply each finite pad's local contact and shoulders.
%% 2. Place the rigid active surfaces at their matching contact planes
cuffFaces_m = zeros(4,3,height(cuffSites)); % Allocate a four-vertex surface for each active pad.
for k = 1:height(cuffSites) % Construct each pad from its measured width and axial position.
    n = normals(k,:); % Read the pad's outward normal.
    tangent = [-n(2),n(1)]; % Read its transverse direction.
    [plane_m,indentation_m] = cuffContactPlane(p,cuffSites(k,:)); % Compute the minimum contact plane from the pad width.
    halfWidth_m = cuffSites.sideXY_m(k)/2; % Read half the finite transverse active width.
    halfLength_m = cuffSites.sideZ_m(k)/2; % Read half the finite axial active length.
    xy = plane_m*n+[-1;1;1;-1]*halfWidth_m*tangent; % Set the face corners on one plane.
    zz = cuffSites.z_m(k)+[-1;-1;1;1]*halfLength_m; % Set the two finite axial pad edges.
    cuffFaces_m(:,:,k) = [xy,zz]; % Store the ideal contact surface without a metal volume.
    cuffIndentation_m(k,1) = indentation_m; % Record the central flattening for this contact.
end % Finish all cuff pads.
nerve.x_m = radius_m .* cos(theta); % Convert transverse radial samples to x.
nerve.y_m = radius_m .* sin(theta); % Convert transverse radial samples to y.
nerve.z_m = repmat(z_m, 1, numel(theta)); % Give every point its axial coordinate.
nerve.cuffFaces_m = cuffFaces_m; % Pass the rigid contacts to the renderer.
nerve.cuffIndentation_m = cuffIndentation_m; % Expose the geometry-derived depths for inspection.
end % Return the nerve surface and finite-area cuff contacts.
