function polys = stationFascicles(q, section) % Fascicle outlines (cell of Nx2, metres) at a station: random shapes after the insertion, or ellipses in legacy mode.
%% Purpose
% One place that tells every figure which fascicle outlines belong to a station. In "insertion" mode these are the random
% irregular shapes displaced by the lead of that station (insertionScenario); otherwise the editable ellipses.
if isfield(q.geometry,'scenario') && ~isempty(q.geometry.scenario) % Insertion mode.
    i = find([q.geometry.scenario.stations.section]==section,1); % Station index.
    if ~isempty(i), polys = q.geometry.scenario.stations(i).polys; return; end % Displaced outlines of this station.
    polys = q.geometry.scenario.polys0; return % Undisturbed outlines.
end
F = q.geometry.fascicles_m; th = linspace(0,2*pi,129)'; th(end) = []; polys = cell(size(F,1),1); % Legacy ellipses.
for k = 1:size(F,1)
    rot = [cosd(F(k,5)) sind(F(k,5));-sind(F(k,5)) cosd(F(k,5))]; % Rotation of the ellipse.
    polys{k} = F(k,1:2)+([cos(th) sin(th)].*F(k,3:4))*rot; % Outline in anatomical coordinates.
end
end
