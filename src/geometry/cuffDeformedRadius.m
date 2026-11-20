function radii_m = cuffDeformedRadius(p, cuffSites, z_m, directions) % Indent the nerve only under each finite cuff pad.
%% 1. Start from the uncompressed nerve outline
z_m = z_m(:); % Treat each axial sample as one row.
base_m = nerveBoundaryRadius(p,directions); % Find the circular or oval outline along each ray.
radii_m = repmat(base_m',numel(z_m),1); % Copy the original contour to all axial levels.
%% 2. Bring the nerve surface onto each rigid contact face
for k = 1:height(cuffSites) % Process each cuff pad independently.
    n = [cosd(cuffSites.theta_deg(k)),sind(cuffSites.theta_deg(k))]; % Read the outward pad normal.
    tangent = [-n(2),n(1)]; % Set the transverse direction across the pad width.
    plane_m = cuffContactPlane(p,cuffSites(k,:)); % Fit the pad plane to the measured width and nerve outline.
    axialDistance_m = abs(z_m-cuffSites.z_m(k))-cuffSites.sideZ_m(k)/2; % Measure distance beyond the pad's axial edges.
    axialFraction = min(1,max(0,axialDistance_m/p.electrode.cuffTransition_m)); % Normalize the axial shoulder width.
    axialWeight = 1-axialFraction.^2.*(3-2*axialFraction); % Blend smoothly to the native nerve beyond the pad.
    projection = directions*n'; % Measure each ray's component toward the pad.
    facing = projection>0; % Select the half of the nerve facing this contact.
    cap_m = plane_m./projection(facing)'; % Find the ray intersection with the fixed contact plane.
    lateral_m = abs(cap_m.*(directions(facing,:)*tangent')'); % Measure position across that plane.
    lateralDistance_m = lateral_m-cuffSites.sideXY_m(k)/2; % Measure distance outside the finite pad width.
    lateralFraction = min(1,max(0,lateralDistance_m/p.electrode.cuffLateralTransition_m)); % Normalize the side shoulder width.
    lateralWeight = 1-lateralFraction.^2.*(3-2*lateralFraction); % Blend smoothly around both pad edges.
    indentation_m = max(0,radii_m(:,facing)-cap_m); % Measure the indentation needed to meet the pad.
    radii_m(:,facing) = radii_m(:,facing)-indentation_m.*(axialWeight*lateralWeight); % Deform only near the active surface.
end % Finish local nerve-to-cuff apposition.
end % Return the surface radius at each requested station and direction.
