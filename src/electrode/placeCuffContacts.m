function [sites, fit] = placeCuffContacts(sites, patient, edgeIndent_m)
%PLACECUFFCONTACTS Place the rigid, flat cuff pads on the (elliptical) nerve surface (D-11).
%   [sites, fit] = placeCuffContacts(sites, patient, edgeIndent_m)
%
%   The cuff pads are rigid and flat: the electrode does not bend, the nerve conforms to it. For each
%   cuff pad:
%     1. landmark: the point of the native nerve outline on the ray at the pad's angle theta_deg
%        (sites.mat gives the angle of every cuff pad; on a circle this is where its centre sits);
%     2. orientation: the face is tangent to the outline at the landmark, i.e. its normal is the
%        outward normal of the ellipse there (on an oval this differs from the radial direction
%        unless the pad sits on one of the two axes);
%     3. depth: the face is pushed along that normal until BOTH side edges of its 200 um width lie
%        edgeIndent_m below the native surface (measured along the normal). The centre of the face
%        is then deeper by the sagitta of the outline under the face. In the 3D geometry the pad box
%        replaces the tissue in front of its face (Boolean union, buildNerveGeometry3D), so the nerve
%        surface under the pad is the flat face itself;
%     4. the box (thickness along the normal) lies outside its face, in the cuff wall.
%   On a circle of radius R this gives back the old placement: centre on the radial line, inner face
%   at R - (edgeIndent + (w/2)^2/(2R)).
%   Shank (intraneural) contacts are not touched: they stay where the tutor's layout puts them.
%
%   Inputs
%     sites         table from loadActiveSites (thickness_m already set to the model value)
%     patient       struct from makePatient (nerveSemiAxes_m, nerveRotation_deg)
%     edgeIndent_m  depth of the face below the native surface at the side edges of the face [m]
%   Outputs
%     sites   same table; for the cuff rows centre_m, r_m, n, w, phi_deg are updated (n points INTO
%             the nerve, as in sites.mat; theta_deg keeps the landmark angle)
%     fit     table, one row per cuff pad: id, theta_deg, centreIndent_m (depth of the face centre
%             below the native surface), edgeIndent_m (min and max over the two edges),
%             tiltFromRadial_deg (angle between the face normal and the radial direction)

arguments
    sites table
    patient (1,1) struct
    edgeIndent_m (1,1) double {mustBePositive}
end

ax = patient.nerveSemiAxes_m;   rot = patient.nerveRotation_deg;
u1 = [cosd(rot) sind(rot)];   u2 = [-sind(rot) cosd(rot)];        % axes of the ellipse
Mq = u1' * u1 / ax(1)^2 + u2' * u2 / ax(2)^2;                       % Q(x) = x * Mq * x'
rows = find(sites.type == "cuff");
nC = numel(rows);
fit = table(sites.id(rows), sites.theta_deg(rows), zeros(nC,1), zeros(nC,2), zeros(nC,1), ...
            'VariableNames', {'id', 'theta_deg', 'centreIndent_m', 'edgeIndent_m', 'tiltFromRadial_deg'});
for j = 1:nC
    k = rows(j);
    ray = [cosd(sites.theta_deg(k)), sind(sites.theta_deg(k))];
    anchor = ray / sqrt(ray * Mq * ray');                           % landmark on the native outline
    n = anchor * Mq;   n = n / norm(n);                             % outward normal there
    t = [-n(2), n(1)];                                              % tangent = width direction
    offset = anchor * t';                                           % tangential coordinate of the landmark
    % native outline along the normal at the two side edges: Q(e n + s t) = 1, outer root e
    A = n * Mq * n';   B = n * Mq * t';   C = t * Mq * t';
    s = offset + [-1 1] * sites.sideXY_m(k) / 2;
    disc = B^2 * s.^2 - A * (C * s.^2 - 1);
    assert(all(disc > 0), 'placeCuffContacts:width', 'Cuff pad %d is wider than the nerve.', sites.id(k));
    eEdge = (-B * s + sqrt(disc)) / A;
    h = min(eEdge) - edgeIndent_m;                                  % plane of the face: x . n = h
    faceCentre = h * n + offset * t;
    centre = faceCentre + sites.thickness_m(k) / 2 * n;             % the box lies outside its face
    nIn = -n;   phi = atan2d(nIn(2), nIn(1));
    sites.centre_m(k, 1:2) = centre;
    sites.r_m(k) = norm(centre);
    sites.phi_deg(k) = phi;
    sites.n(k, :) = [nIn 0];
    sites.w(k, :) = [-sind(phi), cosd(phi), 0];
    fit.centreIndent_m(j) = anchor * n' - h;
    fit.edgeIndent_m(j, :) = sort(eEdge - h);
    fit.tiltFromRadial_deg(j) = acosd(min(1, abs(n * ray')));
end
end
