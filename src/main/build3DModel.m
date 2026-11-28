function M = build3DModel(p, opts)
%BUILD3DMODEL Patient, contacts, geometry, mesh and FE system of the 3D model (no solve).
%   M = build3DModel(p)
%   M = build3DModel(p, opts)
%
%   Purpose
%     Sections 1-4 (up to the assembled system) of main_pipeline_3d, as a function, so that the same
%     model can be built with other parameters (other patients, toy experiments, the mesh-convergence
%     study scripts/run_convergence_3d.m) without touching the production run.
%
%   Two meshers (p.mesh3d.mesher), one solver (our P1 assembly + perineurium thin layer, D-06):
%     "gmsh"        (default since 2026-11-26, D-14) Marco's Gmsh model (src/mesh/meshPersonalDevice.m
%                   + mesh_personal.py): irregular fascicles, one random lead insertion per station
%                   (bent PDMS lead body, displaced fascicles, pads placed along the lead), oval PDMS
%                   sleeve, 5 um elements on the pad edges. Anatomy: src/mesh/gmshAnatomy.m. Needs
%                   Python with the gmsh package (pip install gmsh).
%     "pdetoolbox"  the previous model (buildNerveGeometry3D + generateMesh): straight shanks, circular
%                   fascicles, void cuff wall. Kept as a reference.
%
%   Mesh cache (p.mesh3d.cacheMode). Gmsh does not give the same mesh twice (multithreaded meshing),
%   and every new mesh invalidates the solved contacts. So the mesh is saved once in
%   results/mesh/<patient id>_<mesher>_<key>.mat, where <key> is a fingerprint of everything that
%   defines it (anatomy, lead insertions, contact poses, domain, element sizes; meshKey below), and:
%     "reuse"       (default) load the saved mesh if one exists for the same key, otherwise generate
%                   and save it. Changing the patient, the fascicle layout (toy experiments), the
%                   insertion parameters or the element sizes changes the key: a new mesh is generated
%                   automatically, the old ones stay on disk.
%     "regenerate"  always run the mesher again (to look for differences or errors in the generator);
%                   the new mesh replaces the saved one, and the previous file is kept as
%                   <name>_previous.mat to compare the two.
%
%   Inputs
%     p     parameters (config/pw6_3d_params.m)
%     opts  optional struct:
%       .savePatient  (true)  write data/patients/<id>.mat (the production run does; studies do not)
%       .sizeScale    (1)     multiply every element-size target by this factor (hMin and the edge size
%                             are kept): mesh-convergence study
%       .geom         ([])    pdetoolbox only: a geometry from an earlier call with the same p
%       .verbose      (true)  print the patient table and the mesh statistics
%
%   Output: struct M with the fields used by solve3DContacts and by the plots of main_pipeline_3d:
%     p, patient, R, fasc, sites, cuffFit, cuffAudit, where, geom, nN, nE, elemCell, elemTissue,
%     elemSite, tl (thin-layer mesh: nodes_m, tets, dupOf, ...), nAll, K0, vol, gnd, free, sIdle,
%     root, meshFile

arguments
    p (1,1) struct
    opts (1,1) struct = struct()
end
if ~isfield(opts, 'savePatient'), opts.savePatient = true; end
if ~isfield(opts, 'sizeScale'),   opts.sizeScale = 1; end
if ~isfield(opts, 'geom'),        opts.geom = []; end
if ~isfield(opts, 'verbose'),     opts.verbose = true; end
mm = 1e-3;   um = 1e-6;
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));   % project root
mesher = "pdetoolbox";   if isfield(p.mesh3d, 'mesher'), mesher = string(p.mesh3d.mesher); end
cacheMode = "reuse";     if isfield(p.mesh3d, 'cacheMode'), cacheMode = string(p.mesh3d.cacheMode); end
sizeNames = ["hContact_m" "hEndo_m" "hNerve_m" "hCoarse_m" "hMax_m"];
if opts.sizeScale ~= 1 && mesher == "gmsh"                  % convergence study: scale the size targets
    for n = sizeNames, p.mesh3d.(n) = p.mesh3d.(n) * opts.sizeScale; end   % (pdetoolbox: scaled at meshing)
end

%% ===== 1. PATIENT AND CONTACTS =========================================================
cuffFit = struct();   q = [];
switch mesher
    case "gmsh"
        [q, patient] = gmshAnatomy(p);                       % anatomy + lead insertions of this patient
        [pM, patientM, sites] = personalDevice3D(p, q);      % contacts placed on the bent leads and the oval
        p.personalGeometry = q;                              % kept with the solution (out.params)
        p.electrode.cuffIndent_m = pM.electrode.cuffIndent_m;
        p.domain.nerveLength_m = pM.domain.nerveLength_m;    % = ours (gmshAnatomy sets q.domain)
    case "pdetoolbox"
        patient = makePatient(p.patient.id, p.patient.seed, p);
        sites = loadActiveSites(p.electrode.sitesFile, patient.nerveRadius_m);
        sites.thickness_m(:) = p.electrode.thicknessOverride_m;   % 50 um (slides), not the 5 um of sites.mat
        % rigid flat cuff pads: tangent to the oval, both side edges cuffEdgeIndent_m inside the nerve (D-11)
        [sites, cuffFit] = placeCuffContacts(sites, patient, p.electrode.cuffEdgeIndent_m);
    otherwise
        error('build3DModel:mesher', 'Unknown p.mesh3d.mesher "%s" (gmsh | pdetoolbox).', mesher);
end
R = patient.nerveRadius_m;   fasc = patient.fascicles;   nF = height(fasc);
if opts.savePatient
    genParams = struct('geometry', p.geometry, 'perineurium', p.perineurium, 'mesher', mesher, ...
                       'sitesFile', p.electrode.sitesFile, 'created', string(datetime("now")));
    if mesher == "gmsh", genParams.personalGeometry = q; end
    save(fullfile(root, 'data', 'patients', patient.id + ".mat"), 'patient', 'genParams');
end
if opts.verbose
    fprintf(['Patient %s (seed %d, %s): %s %.3f x %.3f mm (semi-axes, major axis at %g deg; equivalent ' ...
             'R = %.3f mm), %d fascicles; mesher %s\n'], patient.id, patient.seed, patient.level, ...
            patient.nerveShape, patient.nerveSemiAxes_m/mm, patient.nerveRotation_deg, R/mm, nF, mesher);
    disp(table((1:nF)', round(fasc.x_m/mm, 3), round(fasc.y_m/mm, 3), round(fasc.rEndo_m/mm, 3), ...
         round(pi*fasc.rEndo_m.^2/mm^2, 3), 'VariableNames', {'F', 'x_mm', 'y_mm', 'rEq_mm', 'area_mm2'}));
    if mesher == "gmsh"
        sc = q.geometry.scenario;
        for s = sc.stations
            fprintf('   station %+d: lead pierces %s, displaces %s, slides around %s\n', s.section, ...
                    listF(s.pierce), listF(s.displace), listF(s.deflect));
        end
    end
end
where = contactLocation(sites, patient, q);

%% ===== 2-3. GEOMETRY AND MESH (cached) =====================================================
key = meshKey(mesher, p, patient, sites, q, opts.sizeScale);
meshDir = fullfile(root, 'results', 'mesh');   if ~exist(meshDir, 'dir'), mkdir(meshDir); end
meshFile = fullfile(meshDir, sprintf('%s_%s_%s.mat', patient.id, mesher, key));
if cacheMode == "reuse" && isfile(meshFile)
    S = load(meshFile, 'meshData');   md = S.meshData;
    if opts.verbose, fprintf('\nMesh loaded from %s (generated %s)\n', meshFile, md.created); end
else
    if cacheMode == "regenerate" && isfile(meshFile)
        movefile(meshFile, strrep(meshFile, '.mat', '_previous.mat'), 'f');
    end
    md = generateMeshData(mesher, p, patient, sites, q, opts);
    md.created = string(datetime("now"));   md.key = key;
    meshData = rmfieldIf(md, 'g');                           %#ok<NASGU> the CAD geometry is not saved
    save(meshFile, 'meshData', '-v7.3');
    if opts.verbose, fprintf('Mesh saved to %s\n', meshFile); end
end
X = md.X;   T = md.T;   elemCell = md.elemCell;   geom = md.geom;
if isfield(md, 'g'), geom.g = md.g; end
nN = size(X, 2);   nE = size(T, 2);
elemTissue = geom.cellTissue(elemCell);   elemTissue = elemTissue(:)';
elemSite = geom.cellSite(elemCell);       elemSite = elemSite(:)';
cuffAudit = [];   if isfield(md, 'cuffAudit'), cuffAudit = md.cuffAudit; end
if opts.verbose
    fprintf('Mesh: %d nodes, %d linear tetrahedra; per tissue: endo %d, epi %d, muscle %d, contacts %d, PDMS lead %d, PDMS sleeve %d\n', ...
            nN, nE, nnz(elemTissue == 1), nnz(elemTissue == 3), nnz(elemTissue == 4), nnz(elemTissue == 5), ...
            nnz(elemTissue == 6), nnz(elemTissue == 7));
    if ~isempty(cuffAudit), disp(cuffAudit); end
end

%% ===== 4. FE SYSTEM (assembly; the solve is in solve3DContacts) ==========================
sigT = p.materials.sigmaT_Spm;   sigL = p.materials.sigmaL_Spm;   % endo, peri, epi, muscle
Gs = sigT(2) / p.perineurium.realThickness_m;                      % perineurium sheet conductance
if opts.verbose
    fprintf('\nPerineurium thin layer: sheet resistance %.0f Ohm*cm^2 (Gs = %.1f S/m^2)\n', 1e4/Gs, Gs);
end
% thin layer: duplicate the nodes of the fascicle boundaries (not those shared with a contact: there
% the metal shorts the two sides)
contactNodes = unique(T(:, elemTissue == 5));
ifaceNodes = setdiff(md.ifaceNodes, contactNodes);
tl = thinLayerInterface(X, T, elemTissue == 1, ifaceNodes(:), Gs);
nAll = tl.nNodes;
if opts.verbose
    fprintf('   %d interface nodes duplicated, interface area %.1f mm^2 (fascicle perimeters x length: %.1f mm^2)\n', ...
            numel(ifaceNodes), sum(tl.area_m2)/mm^2, sum(2*pi*fasc.rEndo_m) * p.domain.nerveLength_m / mm^2);
end

% conductivity of every element [sigma_x; sigma_y; sigma_z]; codes: 1 endoneurium, 2 (perineurium:
% thin layer, no elements), 3 epineurium, 4 muscle, 5 contacts, 6 PDMS lead body, 7 PDMS sleeve
switch p.electrode.idleContacts
    case "insulating", sIdle = p.electrode.insulatorSigma_Spm;
    case "floating",   sIdle = p.electrode.contactSigma_Spm;
    otherwise, error('Unknown idleContacts "%s".', p.electrode.idleContacts);
end
sP = p.electrode.pdmsSigma_Spm;
tissueSigma = [sigT(1) sigT(1) sigL(1); 0 0 0; sigT(3) sigT(3) sigL(3); sigT(4) sigT(4) sigL(4); ...
               sIdle sIdle sIdle; sP sP sP; sP sP sP]';
sig = tissueSigma(:, elemTissue);
[K0, vol] = assembleConductionP1(tl.nodes_m, tl.tets, sig, p.solve3d.assemblyBlockSize);
K0 = K0 + tl.K;

gnd = md.gndNodes(:);
gnd = unique([gnd; nonzeros(tl.copyOf(gnd))]);
free = true(nAll, 1);   free(gnd) = false;

M = struct('p', p, 'patient', patient, 'R', R, 'fasc', fasc, 'sites', sites, 'cuffFit', cuffFit, ...
           'cuffAudit', cuffAudit, 'where', where, 'geom', geom, 'nN', nN, 'nE', nE, 'elemCell', elemCell, ...
           'elemTissue', elemTissue, 'elemSite', elemSite, 'tl', tl, 'nAll', nAll, 'K0', K0, ...
           'vol', vol, 'gnd', gnd, 'free', free, 'sIdle', sIdle, 'root', root, 'meshFile', meshFile);
end


%% ===== LOCAL FUNCTIONS =================================================================
function md = generateMeshData(mesher, p, patient, sites, q, opts)
%GENERATEMESHDATA Run the mesher; return nodes, tetrahedra, cell roles, interface and ground nodes.
t0 = tic;
switch mesher
    case "gmsh"
        patientM = gmshPatientForMesher(patient, q);
        [geom, X, T, elemCell] = meshPersonalDevice(p, patientM, sites);
        elemTissue = geom.cellTissue(elemCell);   elemTissue = elemTissue(:)';
        elemSite = geom.cellSite(elemCell);       elemSite = elemSite(:)';
        % perineurium = the surface shared by endoneurium and epineurium (where the lead pierces a
        % fascicle the endoneurium touches the PDMS lead instead: no perineurium there)
        ifaceNodes = intersect(unique(T(:, elemTissue == 1)), unique(T(:, elemTissue == 3)));
        % grounded: the outer muscle cylinder, and the end planes if p.domain.endBC = "ground"
        Rout = max(q.geometry.nerveSemiAxes_m) * p.domain.outerRadiusFactor;
        tol = 1e-9;
        grounded = abs(hypot(X(1,:), X(2,:)) - Rout) < tol * max(1, Rout / 1e-3);
        if p.domain.endBC == "ground"
            grounded = grounded | abs(abs(X(3,:)) - p.domain.nerveLength_m / 2) < tol;
        end
        gndNodes = find(grounded)';
        assert(~isempty(gndNodes), 'build3DModel:ground', 'No grounded node found on the outer boundary.');
        cuffAudit = checkCuffMeshContact(X, T, elemTissue, elemSite, sites, q);
        md = struct('X', X, 'T', T, 'elemCell', elemCell, 'geom', geom, 'ifaceNodes', ifaceNodes(:), ...
                    'gndNodes', gndNodes(:), 'cuffAudit', cuffAudit);
    case "pdetoolbox"
        if isempty(opts.geom) || ~isfield(opts.geom, 'g')
            geom = buildNerveGeometry3D(patient, sites, p);
        else
            geom = opts.geom;
        end
        g = geom.g;
        h = geom.faceSize_m * opts.sizeScale;   hs = unique(h(isfinite(h)))';
        Hf = {};
        for s = hs, Hf(end+1:end+2) = {find(h == s)', s}; end %#ok<AGROW>
        meshArgs = {'GeometricOrder', "linear", 'Hmax', p.mesh3d.hMax_m * opts.sizeScale, 'Hmin', p.mesh3d.hMin_m, ...
                    'Hgrad', p.mesh3d.hGrad, 'Hface', Hf};
        if isfield(p.mesh3d, 'hEdge_m') && p.mesh3d.hEdge_m > 0   % refine the contact edges
            contactFaces = idList(arrayfun(@(c) cellFaces(g, c), find(geom.cellTissue == 5), 'UniformOutput', false));
            contactEdges = idList(arrayfun(@(f) faceEdges(g, f), contactFaces, 'UniformOutput', false));
            meshArgs = [meshArgs, {'Hedge', {contactEdges, p.mesh3d.hEdge_m}}];
        end
        gm = generateMesh(g, meshArgs{:});
        msh = gm.Mesh;   X = msh.Nodes;   T = double(msh.Elements);
        elemCell = zeros(1, size(T, 2));
        for c = 1:g.NumCells, elemCell(findElements(msh, "region", "Cell", c)) = c; end
        assert(all(elemCell > 0), 'Some elements belong to no cell.');
        ifaceNodes = findNodes(msh, "region", "Face", geom.interfaceFaces);
        gndNodes = findNodes(msh, "region", "Face", geom.groundFaces);
        md = struct('X', X, 'T', T, 'elemCell', elemCell, 'geom', rmfieldIf(geom, 'g'), ...
                    'ifaceNodes', ifaceNodes(:), 'gndNodes', gndNodes(:), 'g', g);
end
if opts.verbose, fprintf('Mesh generated by %s in %.0f s\n', mesher, toc(t0)); end
end


function pm = gmshPatientForMesher(patient, q)
%GMSHPATIENTFORMESHER The patient struct meshPersonalDevice expects (radius = major semi-axis, fascicles
% with a_m, b_m, angle_deg). With an insertion scenario the mesher uses the outlines of q instead.
F = patient.fascicles;
pm = struct('id', patient.id, 'nerveRadius_m', max(q.geometry.nerveSemiAxes_m), ...
            'fascicles', table(F.x_m, F.y_m, F.rEndo_m, F.rEndo_m, F.rEndo_m, zeros(height(F), 1), ...
                               'VariableNames', {'x_m', 'y_m', 'rEndo_m', 'a_m', 'b_m', 'angle_deg'}));
end


function key = meshKey(mesher, p, patient, sites, q, sizeScale)
%MESHKEY Fingerprint (12 hex digits) of everything that defines the mesh.
m = p.mesh3d;
for f = ["cacheMode"], if isfield(m, f), m = rmfield(m, f); end, end
k = struct('mesher', mesher, 'mesh3d', m, 'sizeScale', sizeScale, 'domain', p.domain, 'cuff', p.cuff, ...
           'contacts', [sites.id, sites.centre_m, sites.n, sites.w, sites.thickness_m, sites.sideXY_m, sites.sideZ_m], ...
           'contactType', string(sites.type)', 'semiAxes', patient.nerveSemiAxes_m, ...
           'rotation', patient.nerveRotation_deg, ...
           'fascicles', [patient.fascicles.x_m, patient.fascicles.y_m, patient.fascicles.rEndo_m]);
if mesher == "gmsh"
    sc = q.geometry.scenario;
    k.outlines = {q.geometry.patient.fascicles.outline};
    k.stations = arrayfun(@(s) struct('path', s.path, 'shift', s.shift, 'z', s.z_m), sc.stations);
    k.leadSizes = [sc.bump_m sc.leadWidth_m sc.leadThickness_m sc.sleeveThickness_m];
    k.cuffOptions = [q.electrode.cuffExtraIndentation_m q.electrode.cuffTransition_m q.electrode.cuffLateralTransition_m];
else
    k.geometry = p.geometry;   k.cuffEdgeIndent = p.electrode.cuffEdgeIndent_m;
end
txt = jsonencode(k);
md5 = java.security.MessageDigest.getInstance('MD5');
md5.update(unicode2native(txt, 'UTF-8'));
h = typecast(md5.digest(), 'uint8');
key = lower(reshape(dec2hex(h(1:6))', 1, []));
end


function where = contactLocation(sites, patient, q)
%CONTACTLOCATION Where each contact lies: fascicle(s) / epineurium crossed by its mid-line. With the
% Gmsh model the fascicle outlines of the contact's station are used (displaced by the lead).
where = strings(height(sites), 1);
for k = 1:height(sites)
    if sites.type(k) == "cuff", where(k) = "cuff, nerve surface"; continue; end
    s = linspace(-0.5, 0.5, 201)' * sites.sideXY_m(k);
    pts = sites.centre_m(k, 1:2) + s * sites.w(k, 1:2);
    lab = zeros(size(pts, 1), 1);
    if ~isempty(q) && isfield(q.geometry, 'scenario') && ~isempty(q.geometry.scenario)
        P = scenarioFascicles(q.geometry.scenario, sites.z_m(k));
        for i = 1:numel(P), lab(inpolygon(pts(:,1), pts(:,2), P{i}(:,1), P{i}(:,2))) = i; end
    else
        F = patient.fascicles;
        for i = 1:height(F), lab(hypot(pts(:,1) - F.x_m(i), pts(:,2) - F.y_m(i)) < F.rEndo_m(i)) = i; end
    end
    parts = strings(0);
    for v = unique(lab)'
        name = "epineurium";   if v > 0, name = "F" + v; end
        parts(end+1) = sprintf('%s %.0f%%', name, 100 * mean(lab == v)); %#ok<AGROW>
    end
    where(k) = strjoin(parts, ", ");
end
end


function s = listF(flags)
%LISTF "F1, F5" for the true entries of a logical vector, "-" if none.
i = find(flags);
if isempty(i), s = "-"; else, s = strjoin("F" + string(i(:)'), ", "); end
end


function s = rmfieldIf(s, f)
%RMFIELDIF Remove field f from struct s if it exists (rmfield errors on a missing field).
if isfield(s, f), s = rmfield(s, f); end
end


function ids = idList(c)
%IDLIST Unique row vector of the ids in a cell array (each cell: a vector, or a cell of vectors).
c = c(:)';
for k = 1:numel(c), if iscell(c{k}), c{k} = [c{k}{:}]; end, c{k} = c{k}(:)'; end
ids = unique([c{:}]);
end
