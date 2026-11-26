function patient = virtualPatient(q, seed, scenario, functionSeed) % One reproducible virtual patient: nerve section, irregular fascicles, functions, motor patches.
%% Purpose
% Virtual patients replace measured data until real sections exist. Seed 1 is the reference patient (V001) used so far: the
% circular layout of config/quicklook3DParams.m with random irregular outlines. Other seeds draw a new layout (repulsion packing
% with the same constraints: 10 fascicles, 5.01 mm2 in total, 1.49 mm2 the largest) and new shapes. Every virtual patient has
% the same structure as a measured one (see docs/patient_data.md), so the rest of the model does not know the difference.
if nargin<2, seed = 1; end, if nargin<3, scenario = "none"; end, if nargin<4, functionSeed = seed; end
ip = insertionParams(); g = q.geometry; q2 = q;
if seed~=1 % New layout for another virtual patient.
    q2.geometry.fascicles_m = packFascicles(q2,RandStream('mt19937ar','Seed',20261200+seed)); % Random layout with the aggregate area constraints (repulsion packing).
end
shapes = randomFascicleShapes(q2,ip,ip.seed+(seed~=1)*seed); % Irregular outlines (roundness from the histology statistics).
n = numel(shapes.polygons); fas = repmat(struct('id',0,'outline',[],'area_m2',0,'roundness',0,'composition',[],'motorFraction',0,'patches',{{}}),n,1);
for k = 1:n, fas(k).id = k; fas(k).outline = shapes.polygons{k}; fas(k).area_m2 = shapes.area_m2(k); fas(k).roundness = shapes.roundness(k); end
patient = struct('id',sprintf('V%03d',seed),'source','virtual','seed',seed,'units','m', ...
    'provenance','Virtual patient: aggregate anatomy from Pusnik 2025 (area 9.20 mm2, 10 fascicles, 5.01 mm2, largest 1.49 mm2), roundness 0.68 +/- 0.08; positions, shapes and functions are random.', ...
    'orientation',struct('x','lateral (assumed)','y','dorsal (assumed)','status','assumed'), ...
    'nerve',struct('outline',shapes.nerve,'semiAxes_m',g.nerveSemiAxes_m,'rotation_deg',g.nerveRotation_deg),'tPeri_m',g.tPeri_m,'fascicles',fas);
if scenario~="none" % A stated hypothesis (random, forearmRule, mrAIN). Default "none": the functions are not known, so none are assigned and no motor patches are drawn.
    patient = assignFunctions(patient,scenario,functionSeed);
else
    cat = functionCatalog(); patient.functionScenario = 'none'; patient.functionSource = 'unknown'; patient.functionIds = cellstr(cat.id)';
end
end
