function drawSection(patient, geom, siteIdsToHighlight, sites)
%DRAWSECTION Nerve, fascicles (with labels) and contact footprints, in mm.
th = linspace(0, 2*pi, 300);   R = patient.nerveRadius_m;   f = patient.fascicles;
if isfield(geom,'nerveOutline'), R=geom.nerveOutline; end % Use the polar radii of the personal outline.
plot(1e3*R.*cos(th), 1e3*R.*sin(th), 'k-', LineWidth=1);   hold on
for i = 1:height(f)
    [xx,yy] = fascicleOutline(f,i,th); % Use the same outline as the finite-element geometry.
    plot(1e3*xx,1e3*yy, '-', ...
         Color=[0.25 0.25 0.25], LineWidth=0.8);
    text(1e3*f.x_m(i), 1e3*f.y_m(i), "F" + i, HorizontalAlignment="center", FontWeight="bold", FontSize=8);
end
for k = 1:numel(geom.footprints)
    if ~isempty(siteIdsToHighlight) && ~any(ismember(geom.footprints(k).siteIds,sites.id)), continue; end % Draw only footprints belonging to this solved station.
    c = geom.footprints(k).corners_m([1:4 1], :);
    hl = any(ismember(geom.footprints(k).siteIds, siteIdsToHighlight));
    col = [0.75 0.75 0.75];   if hl, col = [0.85 0 0]; end
    plot(1e3*c(:,1), 1e3*c(:,2), '-', Color=col, LineWidth=1 + hl);
end
for s = siteIdsToHighlight(:)'
    r = find(sites.id == s);   q = 1e3 * sites.centre_m(r, 1:2);
    d = 0.13 * [cosd(sites.phi_deg(r)), sind(sites.phi_deg(r))];
    text(q(1) + d(1), q(2) + d(2), string(s), Color=[0.85 0 0], FontSize=7, HorizontalAlignment="center");
end
axis equal; box on;
end

