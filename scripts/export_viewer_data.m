%% EXPORT_VIEWER_DATA  Export a solved 3D model for the interactive browser viewer (viewer/index.html).
%
% What it does
%   Reads results/main3d/<id>_3d.mat (written by main_pipeline_3d, nothing is recomputed) and writes
%   results/viewer/<id>/:
%     geometry.js   the model in mm: nerve radius and length, fascicles (x, y, r), cuff, all 48
%                   contacts (centre, orientation, size), which contacts are solved, the sampling grid
%     AS<k>.js      one file per solved contact: the field of that contact sampled on a regular 3D grid
%                   inside the nerve, for a 1 uA stimulus:
%                     V    potential [mV]                  (linear 8-bit code)
%                     E    |E| = |grad V|, field [V/m]     (log 8-bit code, eDecades decades)
%                     J    |J| = |sigma grad V|, current density [A/m^2 = uA/mm^2] (log 8-bit code)
%                     A    activating function for a CATHODIC pulse, -d2V/dz2 over one internode
%                          [V/m^2], endoneurium only, depolarisation (> 0) only (log 8-bit code)
%                     V16  potential [mV per uA], 16-bit LINEAR (int16 x scale): the viewer adds these
%                          up to show ANY pattern of contacts with any currents (superposition)
%   The files are JavaScript (not JSON or binary) because a page opened from disk (file://) may load
%   scripts but not fetch files; the voxel arrays are gzip-compressed and stored as base64 strings
%   (the viewer decompresses them with the browser's DecompressionStream).
%
% Geometry (2026-11-26, Gmsh model D-14): geometry.js also carries the fascicle OUTLINES (irregular,
%   undisplaced) and their displacement by the lead near each station, the PDMS lead body of each
%   station (polygon + axial thickness), the PDMS sleeve (as the cuff wall), and a label for every
%   voxel taken from the mesh (0 = no data, 1 = epineurium, 1+k = fascicle k, 200 = contact,
%   250 = PDMS lead): the viewer uses it for the activating-function mask of the patterns.
%
% How V, E, J are sampled
%   Grid: x, y over the bounding box of the nerve ellipse every gridStepXY_m, z over |z| <= zHalfSpan_m
%   every gridStepZ_m; points outside the nerve ellipse get code 0 = "no data". Each grid point is
%   located (locateRegularGrid, among the nerve and contact elements) in the mesh AFTER the
%   thin-layer node duplication (out.mesh), so a point inside a fascicle takes the endoneurium-side
%   value (see interpolateP1). V is interpolated linearly (P1). E and J are constant in each linear
%   tetrahedron, so they are taken from the element that contains the point: this is the exact FE
%   gradient, and it does not produce the fake spikes a finite difference across the perineurium jump
%   would give.
%
% Activating function (A)
%   The quantity that decides whether a fibre fires (phase 5): for a myelinated fibre along z with
%   nodes of Ranvier every Dz, f = (V(z-Dz) - 2 V(z) + V(z+Dz)) / Dz^2. The solution is for a positive
%   (anodic) current; for a cathodic pulse, the usual way to excite, the sign flips, so A = -f and
%   A > 0 means depolarisation (the fibre is pushed towards firing). Dz = p.plot3d.afSpacing_m (1 mm =
%   internode of a 10 um fibre, internode = 100 D), the same as the preview in main_pipeline_3d. It is
%   computed on the grid lines parallel to z (fibres are straight because the model is extruded) and
%   kept only inside the endoneurium, where the fibres are; hyperpolarised points (A <= 0) get the
%   lowest colour. The last Dz at each end of the grid has no value. This is a map of the stimulus,
%   not yet the threshold criterion (that needs the fibre diameters and an activating-function
%   threshold, phase 5).
%
% Run: setup_paths; export_viewer_data      (then open viewer/index.html in Chrome/Edge/Firefox)
% Time: locating the grid in the mesh the first time (cached in results/viewer/<id>/ for later
% exports of the same run; was ~13 min with pointLocation), then a few seconds per contact.

clear;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);  setup_paths;

%% ===== OPTIONS ========================================================================
patientId    = "UA01";
gridStepXY_m = 50e-6;      % grid spacing across the nerve [m] (contacts are 50 um thick)
gridStepZ_m  = 100e-6;     % grid spacing along the nerve [m] (contact side along z = 200 um)
zHalfSpan_m  = 12.5e-3;    % grid covers |z| <= this [m]; the contacts span |z| <= 10.4 mm
eDecades     = 3;          % |E| and |J| are coded over this many decades below their 99.9th percentile

%% ===== LOAD ===========================================================================
mm = 1e-3;   uA = 1e-6;
S = load(fullfile(root, 'results', 'main3d', patientId + "_3d.mat"), 'out');   out = S.out;
X = double(out.mesh.nodes_m);   T = double(out.mesh.tets);
elemCell = double(out.mesh.elemCell);
elemTissue = out.cells.tissue(elemCell);   elemTissue = elemTissue(:)';
elemSite   = out.cells.site(elemCell);     elemSite = elemSite(:)';
p = out.params;   patient = out.patient;   sites = out.sites;   res = out.results;
R = patient.nerveRadius_m;   F = patient.fascicles;
if isfield(patient, 'nerveSemiAxes_m')                  % elliptical nerve (D-10)
    ax = patient.nerveSemiAxes_m;   rotN = patient.nerveRotation_deg;
else                                                    % result saved before 2026-11-22: circle
    ax = [R R];   rotN = 0;
end
sc = [];                                                % lead insertions of the Gmsh model (D-14)
if isfield(p, 'personalGeometry') && isfield(p.personalGeometry.geometry, 'scenario')
    sc = p.personalGeometry.geometry.scenario;
end
outDir = fullfile(root, 'results', 'viewer', patientId);
if ~exist(outDir, 'dir'), mkdir(outDir); end

%% ===== GRID AND LOCATION IN THE MESH ==================================================
% the grid covers the bounding box of the (rotated) ellipse
extX = hypot(ax(1) * cosd(rotN), ax(2) * sind(rotN));   extY = hypot(ax(1) * sind(rotN), ax(2) * cosd(rotN));
nxh = ceil(extX / gridStepXY_m);   xs = (-nxh:nxh) * gridStepXY_m;
nyh = ceil(extY / gridStepXY_m);   ys = (-nyh:nyh) * gridStepXY_m;
nzh = round(zHalfSpan_m / gridStepZ_m);   zs = (-nzh:nzh) * gridStepZ_m;
[GX, GY] = ndgrid(xs, ys);
inNerve2 = reshape(ellipseQ([GX(:) GY(:)], ax, rotN) < 1, size(GX));
inNerve = repmat(inNerve2, 1, 1, numel(zs));           % x varies fastest = layout of a WebGL 3D texture
clear GX GY
fprintf('Grid %d x %d x %d (%.1f M voxels, %.2f M inside the nerve): locating in the mesh...\n', ...
        numel(xs), numel(ys), numel(zs), numel(inNerve)/1e6, nnz(inNerve)/1e6);
% Locating ~1 M points in ~3 M tetrahedra depends only on the mesh and the grid: it is cached and
% reused as long as both are unchanged (same run). locateRegularGrid searches only the nerve and
% contact elements (the grid points are all inside the nerve).
cacheFile = fullfile(outDir, 'gridLocation_cache.mat');
key = [size(X, 2), size(T, 2), sum(X(:)), numel(xs), numel(ys), numel(zs), gridStepXY_m, gridStepZ_m, 6];
tic;
cached = false;
if isfile(cacheFile)
    C = load(cacheFile);
    cached = isequal(C.key, key) && isfield(C, 'created') && C.created == out.created;
end
if cached
    tetId = C.tetId;   bary = C.bary;
    fprintf('   reused the cached location (%s)\n', cacheFile);
else
    nerveEl = find(ismember(elemTissue, [1 3 5 6]));      % nerve, contacts and PDMS lead (Gmsh model)
    [tl_, bary] = locateRegularGrid(X, T(:, nerveEl), {xs, ys, zs}, inNerve);
    tetId = nan(size(tl_));   tetId(~isnan(tl_)) = nerveEl(tl_(~isnan(tl_)));
    created = out.created;
    save(cacheFile, 'key', 'created', 'tetId', 'bary', '-v7.3');
end
clear C
found = ~isnan(tetId);
fprintf('   done in %.0f s; %d points not found in the mesh (left as "no data")\n', toc, nnz(~found));
% label of every voxel, from the element that contains it: exact also for irregular and displaced
% fascicles (0 = no data, 1 = epineurium, 1+k = fascicle k, 200 = contact, 250 = PDMS lead)
inIdx = find(inNerve);
labels = zeros(size(inNerve), 'uint8');
tF = elemTissue(tetId(found));   cF = out.cells.fascicle(elemCell(tetId(found)));   cF = cF(:)';
lab = zeros(1, nnz(found), 'uint8');
lab(tF == 3) = 1;   lab(tF == 1) = uint8(1 + cF(tF == 1));   lab(tF == 5) = 200;   lab(tF == 6) = 250;
labels(inIdx(found)) = lab;
[ut, ~, loc] = unique(tetId(found));        % elements hit by the grid, and which one each point uses
nodesOfPt = T(:, tetId(found))';            % nFound x 4

% gradients of the barycentric coordinates of the hit elements (same formulas as assembleConductionP1)
Tu = T(:, ut);
p1 = X(:, Tu(1,:));   a = X(:, Tu(2,:)) - p1;   b = X(:, Tu(3,:)) - p1;   c = X(:, Tu(4,:)) - p1;
V6 = dot(a, cross(b, c));
g2 = cross(b, c) ./ V6;   g3 = cross(c, a) ./ V6;   g4 = cross(a, b) ./ V6;   g1 = -(g2 + g3 + g4);

% conductivity of the hit elements [sigma_x; sigma_y; sigma_z] (as in main_pipeline_3d)
sigT = p.materials.sigmaT_Spm;   sigL = p.materials.sigmaL_Spm;
switch p.electrode.idleContacts
    case "insulating", sIdle = p.electrode.insulatorSigma_Spm;
    otherwise,         sIdle = p.electrode.contactSigma_Spm;
end
sP = 1e-7;   if isfield(p.electrode, 'pdmsSigma_Spm'), sP = p.electrode.pdmsSigma_Spm; end
tissueSigma = [sigT(1) sigT(1) sigL(1); 0 0 0; sigT(3) sigT(3) sigL(3); sigT(4) sigT(4) sigL(4); ...
               sIdle sIdle sIdle; sP sP sP; sP sP sP]';      % 6, 7 = PDMS lead body and sleeve (Gmsh model)
sigU = tissueSigma(:, elemTissue(ut));
siteU = elemSite(ut);

%% ===== GEOMETRY FILE ==================================================================
toMm = @(v) round(v / mm, 5);
solved = out.solvedIds(:)';
siteList = cell(height(sites), 1);
for k = 1:height(sites)
    j = find(solved == sites.id(k), 1);
    s = struct('id', sites.id(k), 'type', string(sites.type(k)), 'section', sites.section(k), ...
               'centre', toMm(sites.centre_m(k, :)), 'n', round(sites.n(k, :), 6), ...
               'w', round(sites.w(k, :), 6), 'thickness', toMm(sites.thickness_m(k)), ...
               'sideXY', toMm(sites.sideXY_m(k)), 'sideZ', toMm(sites.sideZ_m(k)), ...
               'solved', ~isempty(j), 'location', "", 'R_kOhm', NaN, 'indent_um', NaN);
    if ~isempty(j)
        s.location = res.location(j);   s.R_kOhm = round(res.R_kOhm(j), 3);
    end
    if isfield(out, 'cuffFit') && isfield(out.cuffFit, 'id') && any(out.cuffFit.id == sites.id(k))   % depth of the pad face centre
        s.indent_um = round(out.cuffFit.centreIndent_m(out.cuffFit.id == sites.id(k)) * 1e6, 2);
    end
    siteList{k} = s;
end
fascList = cell(height(F), 1);
th64 = (0:63)' * 2 * pi / 64;
for i = 1:height(F)
    if isfield(patient, 'outlines')                     % irregular outline (Gmsh model), 64 vertices
        P = patient.outlines{i};   idx = round(linspace(1, size(P, 1) + 1, 65));   P = P(idx(1:end-1), :);
    else                                                % circle
        P = [F.x_m(i) + F.rEndo_m(i) * cos(th64), F.y_m(i) + F.rEndo_m(i) * sin(th64)];
    end
    fascList{i} = struct('id', i, 'x', toMm(F.x_m(i)), 'y', toMm(F.y_m(i)), 'r', toMm(F.rEndo_m(i)), ...
                         'tPeri_um', round(F.tPeri_m(i) * 1e6, 2), 'outline', toMm(P));
end
leadList = {};   disp_ = struct('bump', 0, 'stations', {{}});   cuffKind = "void wall (insulating)";
cuffWall = p.cuff.wallThickness_m;
if ~isempty(sc)                                         % Gmsh model: leads, displacement, PDMS sleeve
    for i = 1:numel(sc.stations)
        st = sc.stations(i);
        pb = polybuffer(st.path, 'lines', sc.leadWidth_m / 2, 'JointType', 'miter');
        pb = intersect(pb, polyshape(sc.nerve));   rg = regions(pb);
        [~, big] = max(arrayfun(@area, rg));   P = rg(big).Vertices;   P = P(all(isfinite(P), 2), :);
        if size(P, 1) > 256, P = P(round(linspace(1, size(P, 1), 256)), :); end   % display only
        leadList{end+1} = struct('section', st.section, 'z', toMm(st.z_m), ...
                                 'thickness', toMm(sc.leadThickness_m), 'xy', toMm(P), ...
                                 'pierced', find(st.pierce(:))', 'displaced', find(st.displace(:))', ...
                                 'slid', find(st.deflect(:))'); %#ok<SAGROW>
    end
    disp_ = struct('bump', toMm(sc.bump_m), 'stations', ...
                   {arrayfun(@(t) struct('z', toMm(t.z_m), 'shift', toMm(t.shift)), sc.stations, 'UniformOutput', false)});
    cuffKind = "PDMS sleeve (insulating)";   cuffWall = sc.sleeveThickness_m;
end
% boundary conditions of the solved run (D-04, D-17), shown in the viewer subtitle: the lateral outer
% surface of the model and the two end planes z = +/- L/2 are each "ground" or "insulated"
bc = struct('lateral', "ground", 'ends', "ground");     % runs saved before the options existed
if isfield(p.domain, 'lateralBC'), bc.lateral = string(p.domain.lateralBC); end
if isfield(p.domain, 'endBC'),     bc.ends    = string(p.domain.endBC);     end
model = struct( ...
    'patient', struct('id', patient.id, 'seed', patient.seed, 'level', string(patient.level), ...
                      'nInner', patient.generatorInfo.nInner), ...
    'units', "mm; fields for a 1 uA monopolar stimulus (return = the grounded outer boundary)", ...
    'bc', bc, ...
    'nerveRadius', toMm(R), 'nerveAxes', toMm(ax), 'nerveRotation', rotN, ...
    'nerveLength', toMm(p.domain.nerveLength_m), ...
    'cuff', struct('length', toMm(p.cuff.length_m), 'wall', toMm(cuffWall), 'kind', cuffKind), ...
    'leads', {leadList}, 'displacement', disp_, 'afSpacing', toMm(p.plot3d.afSpacing_m), ...
    'encoding', "gzip+base64", 'labels', gzB64(labels(:)'), ...
    'outerRadius', toMm(p.domain.outerRadiusFactor * max(ax)), ...
    'fascicles', {fascList}, 'sites', {siteList}, 'solvedIds', solved, ...
    'grid', struct('nx', numel(xs), 'ny', numel(ys), 'nz', numel(zs), 'x0', toMm(xs(1)), ...
                   'y0', toMm(ys(1)), 'z0', toMm(zs(1)), 'dx', toMm(gridStepXY_m), ...
                   'dy', toMm(gridStepXY_m), 'dz', toMm(gridStepZ_m)), ...
    'source', sprintf('results/main3d/%s_3d.mat (created %s)', patientId, out.created), ...
    'exported', string(datetime("now")));
writeJs(fullfile(outDir, 'geometry.js'), 'NV.setModel', model);

%% ===== ONE FIELD FILE PER SOLVED CONTACT ==============================================
inIdx = find(inNerve);
% endoneurium mask and fascicle index of every voxel: from the mesh labels (exact, also where the lead
% displaces a fascicle); fasc2 (one section, for the printout) from the outlines
fascOf = double(labels) - 1;   fascOf(labels < 2 | labels >= 200) = 0;
inEndo = fascOf > 0;
[GX2, GY2] = ndgrid(xs, ys);
fasc2 = zeros(size(GX2));
for i = 1:height(F)
    if isfield(patient, 'outlines')     % Gmsh model: irregular outlines (undisplaced; the displacement by
        P = patient.outlines{i};        % the lead near each station is not shown by the viewer yet)
        fasc2(inpolygon(GX2, GY2, P(:,1), P(:,2))) = i;
    else
        fasc2((GX2 - F.x_m(i)).^2 + (GY2 - F.y_m(i)).^2 < F.rEndo_m(i)^2) = i;
    end
end
nF = height(F);
mAF = round(p.plot3d.afSpacing_m / gridStepZ_m);    % internode in grid steps
dzAF_m = mAF * gridStepZ_m;
fprintf('Activating function over Dz = %.2f mm (p.plot3d.afSpacing_m), endoneurium only\n', dzAF_m / mm);
for j = 1:numel(solved)
    id = solved(j);   tic;
    u = double(out.V_per_A(:, j));
    Vpt = sum(bary(found, :) .* u(nodesOfPt), 2) * uA / mm;            % [mV] for 1 uA
    gradU = g1 .* u(Tu(1,:))' + g2 .* u(Tu(2,:))' + g3 .* u(Tu(3,:))' + g4 .* u(Tu(4,:))';  % [V/m per A]
    sigActive = sigU;
    sigActive(:, siteU == id) = p.electrode.contactSigma_Spm;          % this contact is the conductor
    Emag = vecnorm(gradU, 2, 1)' * uA;                                  % [V/m] for 1 uA
    Jmag = vecnorm(sigActive .* gradU, 2, 1)' * uA;                     % [A/m^2] for 1 uA
    [codeV, qV] = codeLinear(Vpt);
    [codeE, qE] = codeLog(Emag(loc), eDecades);
    [codeJ, qJ] = codeLog(Jmag(loc), eDecades);
    % activating function, cathodic: A = -(V(z-Dz) - 2V(z) + V(z+Dz)) / Dz^2, endoneurium only
    Vfull = nan(size(inNerve));
    Vfull(inIdx(found)) = Vpt * mm;                                     % back to [V] for 1 uA
    A = nan(size(inNerve));
    A(:, :, 1+mAF:end-mAF) = -(Vfull(:, :, 1:end-2*mAF) - 2 * Vfull(:, :, 1+mAF:end-mAF) ...
                               + Vfull(:, :, 1+2*mAF:end)) / dzAF_m^2;   % [V/m^2] for 1 uA
    A(~inEndo) = NaN;
    okA = ~isnan(A);
    [codeA, qA] = codeLog(A(okA), eDecades);
    volA = zeros(size(inNerve), 'uint8');   volA(okA) = codeA;
    qV.unit = "mV";   qE.unit = "V/m";   qJ.unit = "uA/mm^2";   qA.unit = "V/m^2";
    qA.spacing_mm = dzAF_m / mm;
    field = struct('id', id, 'quantities', struct('V', qV, 'E', qE, 'J', qJ, 'A', qA), ...
                   'data', struct('V', toVolume(codeV, inNerve, found), ...
                                  'E', toVolume(codeE, inNerve, found), ...
                                  'J', toVolume(codeJ, inNerve, found), ...
                                  'A', gzB64(volA(:)')), ...
                   'encoding', "gzip+base64");
    % 16-bit linear potential for the patterns (superposition in the viewer): int16 = V / scale * 32767
    v16 = zeros(size(inNerve), 'int16');
    scaleV = max(abs(Vpt));                                              % [mV per uA]
    v16(inIdx(found)) = int16(round(Vpt / scaleV * 32767));
    field.V16 = struct('scale_mV', scaleV, 'data', gzB64(typecast(v16(:)', 'uint8')));
    writeJs(fullfile(outDir, sprintf('AS%d.js', id)), 'NV.addField', field);
    % First look at which fascicles this contact drives. Each grid line along z inside a fascicle is
    % taken as one fibre (uniform density), with its maximum AF along z; the fascicle is summarised by
    % the AF reached by 10% of its fibres (90th percentile), the recruitment level of Zelechowski 2020.
    % The single peak is not used: it comes from the few fibres touching the contact edges, where the
    % activating function is nearly singular. Printout only; the real criterion is phase 5.
    afMaxLine = max(A, [], 3, 'omitnan');                               % nx x ny, NaN outside fascicles
    af10 = zeros(1, nF);
    for i = 1:nF
        a = sort(afMaxLine(fasc2 == i & ~isnan(afMaxLine)), 'descend');
        if ~isempty(a), af10(i) = a(max(1, ceil(0.10 * numel(a)))); end
    end
    [pk, ord] = sort(af10, 'descend');
    fprintf('   AS %2d: V %.2f ... %.2f mV, |E| up to %.0f V/m, |J| up to %.1f A/m^2, AF up to %.0f V/m^2 (per uA)  [%.1f s]\n', ...
            id, qV.lo, qV.hi, qE.hi, qJ.hi, qA.hi, toc);
    fprintf('          AF reached by 10%% of the fibres [V/m^2 per uA]: %s\n', ...
            strjoin(compose('F%d %.0f', [ord(1:4); pk(1:4)]'), ', '));
end
fprintf('Done: %s\nOpen viewer/index.html (add ?patient=%s if it is not UA01).\n', outDir, patientId);


%% ===== LOCAL FUNCTIONS ================================================================
function [code, q] = codeLinear(v)
%CODELINEAR 8-bit code 1...255 of v over [min, max] (0 is reserved for "no data").
lo = min(v);   hi = max(v);
code = uint8(1 + round(254 * (v - lo) / max(hi - lo, eps)));
q = struct('scale', "linear", 'lo', lo, 'hi', hi);
end

function [code, q] = codeLog(v, nDec)
%CODELOG 8-bit code 1...255 of log10(v) over nDec decades below the 99.9th percentile of v.
s = sort(v(v > 0));
hi = s(max(1, round(0.999 * numel(s))));   lo = hi / 10^nDec;
t = (log10(max(v, lo)) - log10(lo)) / nDec;
code = uint8(1 + round(254 * min(max(t, 0), 1)));
q = struct('scale', "log", 'lo', lo, 'hi', hi);
end

function b64 = toVolume(code, inNerve, found)
%TOVOLUME Full voxel array (0 outside the nerve or not found) as a base64 string.
vol = zeros(size(inNerve), 'uint8');
inIdx = find(inNerve);
vol(inIdx(found)) = code;
b64 = gzB64(vol(:)');
end

function b64 = gzB64(bytes)
%GZB64 gzip-compress a byte vector (Java, in memory) and encode it as base64.
bos = java.io.ByteArrayOutputStream();
gz = java.util.zip.GZIPOutputStream(bos);
gz.write(typecast(uint8(bytes(:)'), 'int8'));
gz.close();
b64 = matlab.net.base64encode(typecast(bos.toByteArray()', 'uint8'));
end

function writeJs(fileName, fnName, s)
%WRITEJS Write  fnName(<json>);  to a .js file.
fid = fopen(fileName, 'w');
assert(fid > 0, 'Cannot write %s', fileName);
fprintf(fid, '%s(%s);\n', fnName, jsonencode(s));
fclose(fid);
end
