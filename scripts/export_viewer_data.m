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
%   The files are JavaScript (not JSON or binary) because a page opened from disk (file://) may load
%   scripts but not fetch files; the voxel arrays are base64 strings inside them.
%
% How V, E, J are sampled
%   Grid: x, y from -R to +R every gridStepXY_m, z over |z| <= zHalfSpan_m every gridStepZ_m; points
%   outside the nerve (r > R) get code 0 = "no data". Each grid point is located in the mesh AFTER the
%   thin-layer node duplication (out.mesh), so a point inside a fascicle takes the endoneurium-side
%   value (see interpolateP1). V is interpolated linearly (P1). E and J are constant in each linear
%   tetrahedron, so they are taken from the element that contains the point: this is the exact FE
%   gradient, and it does not produce the fake spikes a finite difference across the perineurium jump
%   would give.
%
% Run: setup_paths; export_viewer_data      (then open viewer/index.html in Chrome/Edge/Firefox)
% Time: locating the grid in the mesh ~13 min the first time (cached in results/viewer/<id>/ for
% later exports of the same run), then a few seconds per contact.

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
outDir = fullfile(root, 'results', 'viewer', patientId);
if ~exist(outDir, 'dir'), mkdir(outDir); end

%% ===== GRID AND LOCATION IN THE MESH ==================================================
nxh = ceil(R / gridStepXY_m);   xs = (-nxh:nxh) * gridStepXY_m;   ys = xs;
nzh = round(zHalfSpan_m / gridStepZ_m);   zs = (-nzh:nzh) * gridStepZ_m;
[GX, GY, GZ] = ndgrid(xs, ys, zs);          % x varies fastest = layout of a WebGL 3D texture
inNerve = GX.^2 + GY.^2 < R^2;
P = [GX(inNerve), GY(inNerve), GZ(inNerve)];
clear GX GY GZ
fprintf('Grid %d x %d x %d (%.1f M voxels, %.2f M inside the nerve): locating in the mesh...\n', ...
        numel(xs), numel(ys), numel(zs), numel(inNerve)/1e6, size(P, 1)/1e6);
% Locating ~1 M points in ~3 M tetrahedra takes minutes (13 min on a laptop), and depends only on
% the mesh and the grid: it is cached and reused as long as both are unchanged (same run).
cacheFile = fullfile(outDir, 'gridLocation_cache.mat');
key = [size(X, 2), size(T, 2), sum(X(:)), numel(xs), numel(zs), gridStepXY_m, gridStepZ_m];
tic;
cached = false;
if isfile(cacheFile)
    C = load(cacheFile);
    cached = isequal(C.key, key);
end
if cached
    tetId = C.tetId;   bary = C.bary;
    fprintf('   reused the cached location (%s)\n', cacheFile);
else
    tr = triangulation(T', X');
    [tetId, bary] = pointLocation(tr, P);
    save(cacheFile, 'key', 'tetId', 'bary', '-v7.3');
end
clear C
found = ~isnan(tetId);
fprintf('   done in %.0f s; %d points not found in the mesh (left as "no data")\n', toc, nnz(~found));
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
tissueSigma = [sigT(1) sigT(1) sigL(1); 0 0 0; sigT(3) sigT(3) sigL(3); sigT(4) sigT(4) sigL(4); ...
               sIdle sIdle sIdle]';
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
               'solved', ~isempty(j), 'location', "", 'R_kOhm', NaN);
    if ~isempty(j)
        s.location = res.location(j);   s.R_kOhm = round(res.R_kOhm(j), 3);
    end
    siteList{k} = s;
end
fascList = cell(height(F), 1);
for i = 1:height(F)
    fascList{i} = struct('id', i, 'x', toMm(F.x_m(i)), 'y', toMm(F.y_m(i)), 'r', toMm(F.rEndo_m(i)), ...
                         'tPeri_um', round(F.tPeri_m(i) * 1e6, 2));
end
model = struct( ...
    'patient', struct('id', patient.id, 'seed', patient.seed, 'level', string(patient.level), ...
                      'nInner', patient.generatorInfo.nInner), ...
    'units', "mm; fields for a 1 uA monopolar stimulus (return = grounded outer boundary)", ...
    'nerveRadius', toMm(R), 'nerveLength', toMm(p.domain.nerveLength_m), ...
    'cuff', struct('length', toMm(p.cuff.length_m), 'wall', toMm(p.cuff.wallThickness_m)), ...
    'outerRadius', toMm(p.domain.outerRadiusFactor * R), ...
    'fascicles', {fascList}, 'sites', {siteList}, 'solvedIds', solved, ...
    'grid', struct('nx', numel(xs), 'ny', numel(ys), 'nz', numel(zs), 'x0', toMm(xs(1)), ...
                   'y0', toMm(ys(1)), 'z0', toMm(zs(1)), 'dx', toMm(gridStepXY_m), ...
                   'dy', toMm(gridStepXY_m), 'dz', toMm(gridStepZ_m)), ...
    'source', sprintf('results/main3d/%s_3d.mat (created %s)', patientId, out.created), ...
    'exported', string(datetime("now")));
writeJs(fullfile(outDir, 'geometry.js'), 'NV.setModel', model);

%% ===== ONE FIELD FILE PER SOLVED CONTACT ==============================================
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
    qV.unit = "mV";   qE.unit = "V/m";   qJ.unit = "uA/mm^2";
    field = struct('id', id, 'quantities', struct('V', qV, 'E', qE, 'J', qJ), ...
                   'data', struct('V', toVolume(codeV, inNerve, found), ...
                                  'E', toVolume(codeE, inNerve, found), ...
                                  'J', toVolume(codeJ, inNerve, found)));
    writeJs(fullfile(outDir, sprintf('AS%d.js', id)), 'NV.addField', field);
    fprintf('   AS %2d: V %.2f ... %.2f mV, |E| up to %.0f V/m, |J| up to %.1f A/m^2 (per uA)  [%.1f s]\n', ...
            id, qV.lo, qV.hi, qE.hi, qJ.hi, toc);
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
b64 = matlab.net.base64encode(vol(:)');
end

function writeJs(fileName, fnName, s)
%WRITEJS Write  fnName(<json>);  to a .js file.
fid = fopen(fileName, 'w');
assert(fid > 0, 'Cannot write %s', fileName);
fprintf(fid, '%s(%s);\n', fnName, jsonencode(s));
fclose(fid);
end
