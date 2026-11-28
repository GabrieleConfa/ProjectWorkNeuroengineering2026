function [plane_m, indentation_m, contactWidth_m, measured, frame] = cuffContactPlane(p, site) % Fit one rigid flat pad to the nerve's native outline.
%% Purpose
% A cuff (external) contact is a rigid flat pad pressed on the oval nerve. This function finds the plane of its active face:
% tangent to the ellipse at the pad's polar landmark (site.theta_deg) and pushed in until BOTH side edges of the pad touch
% the nerve, plus cuffExtraIndentation_m (or by a measured depth, if p.electrode.cuffMeasuredDepth_m has a row for this pad).
% The ellipse is x'Qx = 1 with semi-axes p.geometry.nerveSemiAxes_m rotated by p.geometry.nerveRotation_deg; the line/ellipse
% intersections are solved analytically in the pad frame (n = outward normal, t = tangent across the pad).
%% Inputs / outputs
%   p               parameter struct (quicklook3DParams): p.geometry.nerveSemiAxes_m [m], nerveRotation_deg [deg],
%                   p.electrode.cuffExtraIndentation_m [m], optional p.electrode.cuffMeasuredDepth_m ([id depth_m] rows)
%   site            one row of the contact table (loadActiveSites): theta_deg, sideXY_m (pad width) [m], id
%   plane_m         distance of the pad face from the nerve centre along n [m]
%   indentation_m   depth of the face below the undeformed surface at the pad centre [m]
%   contactWidth_m  width of the pad that touches the undeformed nerve [m] (<= sideXY_m)
%   measured        true if the depth came from cuffMeasuredDepth_m, false if from the full-contact assumption
%   frame           struct: n, t (unit vectors), offset (tangential position of the landmark [m]), centre (face centre
%                   [m]), cornerAngles (polar angles of the two pad edges [rad], used by the mesher and the plots)
%% 1. Resolve the rotated oval in the pad's local directions
ray = [cosd(site.theta_deg),sind(site.theta_deg)]; % Keep the supplied polar location as the native contact landmark.
a=p.geometry.nerveSemiAxes_m(1); b=p.geometry.nerveSemiAxes_m(2); % Read the anatomical semi-axes.
c=cosd(p.geometry.nerveRotation_deg); s=sind(p.geometry.nerveRotation_deg); % Rotate the quadratic form into world coordinates.
rotation=[c,s;-s,c]; Q=rotation'*diag([1/a^2 1/b^2])*rotation; % Represent the ellipse by x'Qx = 1.
anchor=ray/sqrt(ray*Q*ray'); % Intersect the nominal polar ray with the native oval.
n=anchor*Q; n=n/norm(n); t=[-n(2),n(1)]; % Use the normalized surface gradient and its perpendicular tangent.
offset=dot(anchor,t); % Retain the tangential coordinate of the landmark during normal closure.
nLocal = rotation*n'; % Express the pad normal in oval coordinates.
tLocal = rotation*t'; % Express the transverse pad direction in oval coordinates.
A = (nLocal(1)/a)^2+(nLocal(2)/b)^2; % Compute the quadratic coefficient along the normal.
B = nLocal(1)*tLocal(1)/a^2+nLocal(2)*tLocal(2)/b^2; % Compute the normal-lateral coupling.
C = (tLocal(1)/a)^2+(tLocal(2)/b)^2; % Compute the quadratic coefficient across the pad.
%% 2. Locate a rigid plane under both transverse pad edges
lateral_m = offset+site.sideXY_m/2*[-1,1]; % Evaluate the two measured active-face edges.
discriminant = B^2*lateral_m.^2+A*(1-C*lateral_m.^2); % Solve the line-oval intersections analytically.
assert(all(discriminant>=0), 'Cuff pad is wider than the nerve at this location.'); % Reject an impossible contact footprint.
edgeNormal_m = (-B*lateral_m+sqrt(discriminant))/A; % Find the exterior nerve point at each lateral edge.
nativeCentre_m = dot(anchor,n); % Find the original surface below the centre of the pad.
plane_m = min(edgeNormal_m)-p.electrode.cuffExtraIndentation_m; % Indent enough for the whole rigid pad to fit.
measured=false; % Use the geometric full-contact assumption until a measurement exists for this pad.
if isfield(p.electrode,'cuffMeasuredDepth_m') && ~isempty(p.electrode.cuffMeasuredDepth_m) % Allow device-specific indentation data.
    data=p.electrode.cuffMeasuredDepth_m; % Read absolute centre indentations in metres.
    assert(size(data,2)==2 && all(isfinite(data),'all') && all(data(:,2)>=0) && numel(unique(data(:,1)))==size(data,1),'Invalid cuff indentation measurements.'); % Reject ambiguous or invalid rows.
    value=data(data(:,1)==site.id,2); % Match the physical electrode ID.
    if ~isempty(value), plane_m=nativeCentre_m-value; measured=true; end % Apply the observed indentation even for partial contact.
end % Finish selecting the plane.
indentation_m = nativeCentre_m-plane_m; % Return the resulting central flattening depth.
assert(plane_m>0 && indentation_m>=0, 'Cuff plane must lie inside the original nerve contour.'); % Reject an invalid fitted plane.
disc=(B*plane_m)^2-C*(A*plane_m^2-1); % Intersect the native ellipse with the rigid face to obtain its contact interval.
roots=(-B*plane_m+[-1 1]*sqrt(max(0,disc)))/C; % Express both intersections along the pad tangent.
contactWidth_m=max(0,min(lateral_m(2),roots(2))-max(lateral_m(1),roots(1))); % Keep only the contacting portion of the finite width.
frame.n=n; frame.t=t; frame.offset=offset; frame.centre=plane_m*n+offset*t; % Locate the active face at its actual point on the oval.
corners=frame.centre+[-1;1]*(site.sideXY_m/2)*t; % Locate both transverse pad edges.
frame.cornerAngles=atan2(corners(:,2),corners(:,1)); % Supply the exact meshing and plotting rays.
end % Return the imposed displacement and geometric footprint.
