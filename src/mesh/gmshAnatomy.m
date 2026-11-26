function [q, patient] = gmshAnatomy(p)
%GMSHANATOMY Anatomy and lead insertions of the Gmsh model, for the patient chosen in p.patient.
%   [q, patient] = gmshAnatomy(p)
%
%   Purpose
%     The Gmsh mesher (src/mesh/meshPersonalDevice.m + mesh_personal.py, from Marco's branch, D-14)
%     reads its anatomy from Marco's parameter struct q (config/quicklook3DParams.m): nerve outline,
%     fascicle outlines (irregular), and one random lead insertion per station (q.geometry.scenario,
%     src/insertion/insertionScenario.m). This function builds that q for the patient of THIS run,
%     and also returns the patient in the format of our pipeline (plots, viewer, block study).
%
%   Patient sources (p.patient.source), so that many patients and toy experiments can be run with
%   the same code (only the source, seed or file changes):
%     "generator"  our upper-arm generator makePatient(p.patient.id, p.patient.seed, p) (D-09,
%                  literature-based layout); fascicle outlines from it: circles, or random irregular
%                  shapes with the same centres and areas (p.patient.fascicleShape = "irregular",
%                  randomFascicleShapes, roundness statistics of Pusnik 2025)
%     "virtual"    Marco's virtual patient V<seed> (virtualPatient: repulsion packing + irregular
%                  shapes; seed 1 = his reference V001)
%     "file"       a patient from a JSON file (p.patient.file; format docs/patient_data.md), e.g. a
%                  measured section
%   The lead insertions are then drawn on that anatomy with the model chosen in
%   config/insertionParams.m (pierceModel), with fixed seeds: reproducible.
%
%   Inputs
%     p   parameters (config/pw6_3d_params.m): p.patient (id, seed, source, fascicleShape, file),
%         p.domain.nerveLength_m, p.perineurium.realThickness_m
%
%   Outputs
%     q        Marco's struct with q.geometry.patient (outlines) and q.geometry.scenario (insertions)
%     patient  our patient struct (id, seed, level, nerveShape, nerveSemiAxes_m, nerveRotation_deg,
%              nerveRadius_m, fascicles table x_m/y_m/rEndo_m (circular equivalents), tPeri_m,
%              generatorInfo) plus .outlines (cell of undisplaced fascicle outlines [m]) and
%              .nerveOutline_m (polygon)

arguments
    p (1,1) struct
end
src = "generator";      if isfield(p.patient, 'source'), src = string(p.patient.source); end
shape = "irregular";    if isfield(p.patient, 'fascicleShape'), shape = string(p.patient.fascicleShape); end

switch src
    case "virtual"
        q = quicklook3DParams(struct('virtualSeed', p.patient.seed));
    case "file"
        q = quicklook3DParams(struct('patientFile', string(p.patient.file)));
    case "generator"
        q = quicklook3DParams();                              % Marco's defaults; anatomy replaced below
        pt = makePatient(p.patient.id, p.patient.seed, p);
        F = pt.fascicles;
        q.geometry.nerveSemiAxes_m = pt.nerveSemiAxes_m;
        q.geometry.nerveRotation_deg = pt.nerveRotation_deg;
        q.geometry.fascicles_m = [F.x_m, F.y_m, F.rEndo_m, F.rEndo_m, zeros(height(F), 1)];
        q.geometry.tPeri_m = p.perineurium.realThickness_m;
        nerve = nerveOutline(pt.nerveSemiAxes_m, pt.nerveRotation_deg, 720);
        if shape == "irregular"                               % same centres and areas, random shapes
            ip = insertionParams();
            shapes = randomFascicleShapes(q, ip, ip.seed + p.patient.seed);
            outl = shapes.polygons;
        else                                                  % the circles of the generator
            th = (0:127)' * 2 * pi / 128;
            outl = arrayfun(@(k) [F.x_m(k) + F.rEndo_m(k) * cos(th), F.y_m(k) + F.rEndo_m(k) * sin(th)], ...
                            (1:height(F))', 'UniformOutput', false);
        end
        mp = patientFromOutlines(p.patient.id, nerve, outl, q.geometry.tPeri_m);
        for k = 1:numel(outl)                                 % roundness 4A/(pi*Feret^2), as in the histology
            P = outl{k};   feret = max(hypot(P(:,1) - P(:,1)', P(:,2) - P(:,2)'), [], 'all');   % no toolbox needed
            mp.fascicles(k).roundness = 4 * mp.fascicles(k).area_m2 / (pi * feret^2);
        end
        mp.source = 'generator';   mp.seed = p.patient.seed;
        mp.nerve.semiAxes_m = pt.nerveSemiAxes_m;   mp.nerve.rotation_deg = pt.nerveRotation_deg;
        validatePatient(mp, q);
        q.geometry.patient = mp;
        q.geometry.scenario = insertionScenario(q);           % lead insertions on THIS anatomy
    otherwise
        error('gmshAnatomy:source', 'Unknown p.patient.source "%s" (generator | virtual | file).', src);
end
q.domain.length_m = p.domain.nerveLength_m;                   % keep our domain length (D-02)

% --- our patient struct (for plots, viewer, block study) -----------------------------------------------
mp = q.geometry.patient;
outl = {mp.fascicles.outline}';
area = cellfun(@(P) polyarea(P(:,1), P(:,2)), outl);
ctr = cell2mat(cellfun(@(P) mean(P, 1), outl, 'UniformOutput', false));
F = table(ctr(:,1), ctr(:,2), sqrt(area / pi), repmat(q.geometry.tPeri_m, numel(outl), 1), ...
          'VariableNames', {'x_m', 'y_m', 'rEndo_m', 'tPeri_m'});
if src == "generator"
    patient = pt;                                             % keep everything the generator gives
    patient.fascicles.x_m = F.x_m;   patient.fascicles.y_m = F.y_m;   patient.fascicles.rEndo_m = F.rEndo_m;
else
    ax = q.geometry.nerveSemiAxes_m;
    patient = struct('id', string(mp.id), 'seed', mp.seed, 'level', "upper arm (" + string(mp.source) + ")", ...
                     'nerveShape', "ellipse", 'nerveSemiAxes_m', ax, 'nerveRotation_deg', q.geometry.nerveRotation_deg, ...
                     'nerveRadius_m', sqrt(prod(ax)), 'fascicles', F, ...
                     'generatorInfo', struct('nInner', NaN, 'targetInner', NaN, 'restarts', 0));
end
patient.outlines = outl;
patient.nerveOutline_m = mp.nerve.outline;
patient.fascicleShape = shape;
patient.source = src;
end
