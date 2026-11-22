function patient = makePatient(id, seed, params)
%MAKEPATIENT Create a virtual patient: nerve radius and fascicles of the cross-section.
%   patient = makePatient(id, seed, params)
%
%   Contract: docs/architecture.md §4.
%
%   Inputs
%     id      string, e.g. "UA01"
%     seed    integer for rng(seed); NaN -> the tutor's fixed section
%             (params.geometry.tutorSection), i.e. patient P00
%     params  parameter struct (default_params, pw6_3d_params, ...). The generated case uses
%             params.geometry.generator, params.perineurium and the shank-contact positions
%             (params.electrode.sitesFile, params.electrode.thicknessOverride_m).
%
%   Output: struct
%     id, seed, level
%     nerveRadius_m   equivalent nerve radius sqrt(area/pi) [m] (the radius for a circle)
%     nerveShape      "ellipse" | "circle" (D-10)
%     nerveSemiAxes_m [a b] [m], a along nerveRotation_deg (a = b = nerveRadius_m for a circle)
%     nerveRotation_deg  angle of the first semi-axis from x [deg]
%     fascicles       table: x_m, y_m, rEndo_m, tPeri_m  [m] (anatomical values; the 3D model may
%                     draw a thicker perineurium, see params.perineurium)
%     generatorInfo   restarts / iterations / inner fascicles (generated case only)
%
%   The patient is regenerated exactly from (seed, params.geometry, params.perineurium, sites.mat).
%   Save with: save(fullfile('data','patients', id + ".mat"), 'patient', 'genParams')

arguments
    id (1,1) string
    seed (1,1) double
    params (1,1) struct
end

if isnan(seed)                                        % tutor's fixed section P00
    ts = params.geometry.tutorSection;
    f = ts.fascicles_m;
    patient = struct('id', id, 'seed', seed, 'level', "tutor section P00", ...
        'nerveRadius_m', ts.nerveRadius_m, 'nerveShape', "circle", ...
        'nerveSemiAxes_m', ts.nerveRadius_m * [1 1], 'nerveRotation_deg', 0, ...
        'fascicles', table(f(:,1), f(:,2), f(:,3), repmat(ts.tPeri_m, size(f,1), 1), ...
                           'VariableNames', {'x_m', 'y_m', 'rEndo_m', 'tPeri_m'}), ...
        'generatorInfo', struct());
    return
end

gen = params.geometry.generator;
R = sqrt(gen.nerveArea_m2 / pi);
sites = loadActiveSites(params.electrode.sitesFile, R);
footprints = shankFootprints(sites, params.electrode.thicknessOverride_m);

rng(seed, 'twister');
[R, fasc, info, nerve] = generateFascicleLayout(gen, params.perineurium, footprints);
patient = struct('id', id, 'seed', seed, 'level', string(gen.level), 'nerveRadius_m', R, ...
                 'nerveShape', nerve.shape, 'nerveSemiAxes_m', nerve.semiAxes_m, ...
                 'nerveRotation_deg', nerve.rotation_deg, 'fascicles', fasc, 'generatorInfo', info);
end

function fp = shankFootprints(sites, thickness_m)
%SHANKFOOTPRINTS Rectangles [m] (4x2 corners) of the distinct shank contacts in the x-y plane.
% Sections +k and -k share the same footprint, so duplicates are removed.
S = sites(sites.type == "intraneural", :);
key = round([S.centre_m(:,1:2), S.phi_deg / 1e3] * 1e7);
[~, iu] = unique(key, 'rows', 'stable');
fp = cell(numel(iu), 1);
sgn = [-1 -1; 1 -1; 1 1; -1 1];
for k = 1:numel(iu)
    s = S(iu(k), :);
    fp{k} = s.centre_m(1:2) + sgn(:,1) * (thickness_m/2) * s.n(1:2) + sgn(:,2) * (s.sideXY_m/2) * s.w(1:2);
end
end
