function polys = scenarioFascicles(sc, z_m) % Fascicle outlines of the scenario at an axial position z (cell of Nx2, metres).
%% Purpose
% At a station the lead displaces some fascicles. The displacement fades linearly to zero over sc.bump_m on both sides
% of the station plane, so the 3D fascicles are continuous. Away from the stations the undisturbed outlines are used.
% Inputs: sc  insertion scenario (insertionScenario): polys0 (undisturbed outlines [m]), stations(:).z_m [m] and .shift
%             (n x 2 displacement of every fascicle [m]), bump_m (axial half-length of the displacement [m]);
%         z_m axial position [m].
% Output: polys, n x 1 cell of N x 2 outlines [m] at z_m.
polys =sc.polys0; % Start from the undisturbed shapes.
for s = sc.stations % Add the displacement of every station.
    w = max(0,1-abs(z_m-s.z_m)/sc.bump_m); % Triangular weight: 1 at the station plane, 0 beyond the bump length.
    if w==0, continue; end
    for k = 1:numel(polys), polys{k} = polys{k}+w*s.shift(k,:); end % Translate fascicle k by the weighted displacement.
end
end
