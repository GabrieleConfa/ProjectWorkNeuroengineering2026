%% RUN_PIPELINE  Full pipeline for one virtual patient (skeleton).
% Owner: WP5. Calls the modules in the order of docs/architecture.md §3.
% Each step stops with "not implemented" until its owner fills it in.
% Run from the project root:   setup_paths; run_pipeline

params = default_params();
sites  = loadActiveSites(params.electrode.sitesFile, ...
                         params.geometry.tutorSection.nerveRadius_m);

% --- Phase 1: patient and cross-section (WP1) ---
patient = makePatient("P00", NaN, params);              % tutor's fixed section
xsec    = buildCrossSection2D(patient, params);

% --- Phase 3: 3D geometry, electrode, mesh (WP2) ---
geom3      = buildNerveGeometry3D(xsec, params);
geom3      = addActiveSites(geom3, sites, params);
nerveModel = meshNerveModel(geom3, patient, sites, params);

% --- Phases 4-5: one solve per AS, then activation (WP3, WP4) ---
fibres = makeFibres(patient, params, 1);
thr = table();
for siteId = 1:height(sites)
    sol = solveUnitCurrent(nerveModel, siteId, params);
    af  = activatingFunction(@(x,y,z) interpolateSolution(sol.result, x, y, z), fibres);
    thr = [thr; thresholdPerFascicle(af, fibres, siteId, params)]; %#ok<AGROW>
end

% --- Safety (WP5) ---
safety = checkSafety(thr, sites, params);

outDir = fullfile('results', 'pipeline');
if ~isfolder(outDir), mkdir(outDir); end
save(fullfile(outDir, patient.id + "_results.mat"), 'patient', 'thr', 'safety', 'params');
