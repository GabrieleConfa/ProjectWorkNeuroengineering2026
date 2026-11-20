%% MAIN_PIPELINE_3D  3D median-nerve model, upper arm: patient -> geometry -> mesh -> FEM solution.
%
% Purpose
%   Generates a virtual patient (10 fascicles, mid-upper arm), builds the 3D model of the nerve with
%   the cuff and the 48 contacts of the electrode, meshes it and solves the volume-conductor problem
%   div(sigma grad V) = 0 for a unit current injected by each selected contact (monopolar, return =
%   grounded outer surface). The problem is linear: the solution for any current I is I times it.
%
% Inputs
%   config/pw6_3d_params.m      ALL parameters and modelling choices (change things there)
%   data/electrode/sites.mat    contact layout provided by the tutor (via loadActiveSites)
%
% Outputs
%   data/patients/<id>.mat            the virtual patient (fascicles) + the parameters that made it
%   results/main3d/<id>_*.mat         mesh, cell roles and potentials (struct "out")
%   results/main3d/<id>_*.png         cross-section, z layout, potential maps and profiles
%
% Model in short (choices and reasons: docs/decisions.md D-02 ... D-09)
%   - geometry: circular nerve (the rigid cuff is tight on it), fascicles generated from the
%     literature (makePatient), the same section along the 70 mm (layered extrusion);
%   - cuff: 50 mm, centred on the contacts, wall = void (perfect insulator); flat cuff contacts
%     pressed 10 um into the nerve; shank contacts only (no shank body yet);
%   - tissues: endoneurium (anisotropic, sigma_z = longitudinal), epineurium, muscle (anisotropic);
%     PERINEURIUM = thin layer on the fascicle boundary (potential jump, sheet resistance rho*t);
%   - active contact: high-conductivity cell with a uniform source I/volume (equipotential, carries
%     exactly I); idle contacts insulating (as in 2D);
%   - outer muscle surface (and the end planes, by default) grounded.
%   The FE system is assembled by the project's own P1 code (src/fem/), checked by
%   tests/testThinLayerFEM.m; the PDE Toolbox builds the geometry and the mesh.
%
% Units: SI inside the code (m, S/m, A, V); mm, um, uA, mV only in prints and plots.
% Run: press Run in the editor, or type main_pipeline_3d (the script adds the project paths itself).
% Time: geometry ~ minutes (many cell merges), mesh ~ minutes, each contact ~ 1 solve.

clear; close all; clc;
mm = 1e-3;   um = 1e-6;
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));   % project root
addpath(root);  setup_paths;
outDir = fullfile(root, 'results', 'main3d');
if ~exist(outDir, 'dir'), mkdir(outDir); end


%% ===== 1. PARAMETERS, PATIENT, CONTACTS ==============================================
p = pw6_3d_params();
patient = makePatient(p.patient.id, p.patient.seed, p);
R = patient.nerveRadius_m;   fasc = patient.fascicles;   nF = height(fasc);
genParams = struct('geometry', p.geometry, 'perineurium', p.perineurium, ...
                   'sitesFile', p.electrode.sitesFile, 'created', string(datetime("now")));
save(fullfile(root, 'data', 'patients', patient.id + ".mat"), 'patient', 'genParams');
fprintf('Patient %s (seed %d, %s): R = %.3f mm, %d fascicles, %d inner (target %d), %d restarts\n', ...
        patient.id, patient.seed, patient.level, R/mm, nF, patient.generatorInfo.nInner, ...
        patient.generatorInfo.targetInner, patient.generatorInfo.restarts);
disp(table((1:nF)', round(fasc.x_m/mm, 3), round(fasc.y_m/mm, 3), round(fasc.rEndo_m/mm, 3), ...
     round(pi*fasc.rEndo_m.^2/mm^2, 3), 'VariableNames', {'F', 'x_mm', 'y_mm', 'rEndo_mm', 'area_mm2'}));

sites = loadActiveSites(p.electrode.sitesFile, R);
sites.thickness_m(:) = p.electrode.thicknessOverride_m;   % 50 um (slides), not the 5 um of sites.mat
isCuff = sites.type == "cuff";
% flat cuff contacts pressed into the nerve: inner face at R - indent
sites.r_m(isCuff) = R - p.electrode.cuffIndent_m + sites.thickness_m(isCuff) / 2;
sites.centre_m(isCuff, 1:2) = sites.r_m(isCuff) .* [cosd(sites.theta_deg(isCuff)), sind(sites.theta_deg(isCuff))];
where = contactLocation(sites, fasc);


%% ===== 2. GEOMETRY ===================================================================
geom = buildNerveGeometry3D(patient, sites, p);
g = geom.g;
fprintf('\nGeometry: %d cells (endo %d, epi %d, muscle %d, contacts %d), %d faces, %d merges, %.0f s\n', ...
        g.NumCells, nnz(geom.cellTissue == 1), nnz(geom.cellTissue == 3), nnz(geom.cellTissue == 4), ...
        nnz(geom.cellTissue == 5), g.NumFaces, geom.nMerges, geom.time_s);


%% ===== 3. MESH =======================================================================
h = geom.faceSize_m;   hs = unique(h(isfinite(h)))';
Hf = {};
for s = hs, Hf(end+1:end+2) = {find(h == s)', s}; end %#ok<AGROW>
tic;
gm = generateMesh(g, GeometricOrder="linear", Hmax=p.mesh3d.hMax_m, Hmin=p.mesh3d.hMin_m, ...
                  Hgrad=p.mesh3d.hGrad, Hface=Hf);
msh = gm.Mesh;   X = msh.Nodes;   T = double(msh.Elements);
nN = size(X, 2);   nE = size(T, 2);
elemCell = zeros(1, nE);
for c = 1:g.NumCells, elemCell(findElements(msh, "region", "Cell", c)) = c; end
assert(all(elemCell > 0), 'Some elements belong to no cell.');
elemTissue = geom.cellTissue(elemCell)';   elemSite = geom.cellSite(elemCell)';
fprintf('Mesh: %d nodes, %d linear tetrahedra (%.0f s); per tissue: endo %d, epi %d, muscle %d, contacts %d\n', ...
        nN, nE, toc, nnz(elemTissue == 1), nnz(elemTissue == 3), nnz(elemTissue == 4), nnz(elemTissue == 5));


%% ===== 4. FE SYSTEM AND SOLUTION =====================================================
sigT = p.materials.sigmaT_Spm;   sigL = p.materials.sigmaL_Spm;   % endo, peri, epi, muscle
Gs = sigT(2) / p.perineurium.realThickness_m;                      % perineurium sheet conductance
fprintf('\nPerineurium thin layer: sheet resistance %.0f Ohm*cm^2 (Gs = %.1f S/m^2)\n', 1e4/Gs, Gs);

% thin layer: duplicate the nodes of the fascicle boundaries (not those shared with a contact)
contactNodes = unique(T(:, elemTissue == 5));
ifaceNodes = setdiff(findNodes(msh, "region", "Face", geom.interfaceFaces), contactNodes);
tl = thinLayerInterface(X, T, elemTissue == 1, ifaceNodes(:), Gs);
nAll = tl.nNodes;
fprintf('   %d interface nodes duplicated, interface area %.1f mm^2 (fascicle perimeters x length: %.1f mm^2)\n', ...
        numel(ifaceNodes), sum(tl.area_m2)/mm^2, sum(2*pi*fasc.rEndo_m) * p.domain.nerveLength_m / mm^2);

% conductivity of every element [sigma_x; sigma_y; sigma_z]
switch p.electrode.idleContacts
    case "insulating", sIdle = p.electrode.insulatorSigma_Spm;
    case "floating",   sIdle = p.electrode.contactSigma_Spm;
    otherwise, error('Unknown idleContacts "%s".', p.electrode.idleContacts);
end
tissueSigma = [sigT(1) sigT(1) sigL(1); 0 0 0; sigT(3) sigT(3) sigL(3); sigT(4) sigT(4) sigL(4); sIdle sIdle sIdle]';
sig = tissueSigma(:, elemTissue);
[K0, vol] = assembleConductionP1(tl.nodes_m, tl.tets, sig);
K0 = K0 + tl.K;

gnd = findNodes(msh, "region", "Face", geom.groundFaces);
gnd = unique([gnd(:); nonzeros(tl.copyOf(gnd))]);
free = true(nAll, 1);   free(gnd) = false;

if isempty(p.solve3d.sections), solveIds = sites.id;
else, solveIds = sites.id(ismember(sites.section, p.solve3d.sections));
end
nS = numel(solveIds);   I0 = p.stim.unitCurrent_A;
if p.solve3d.saveSingle, V = zeros(nAll, nS, 'single'); else, V = zeros(nAll, nS); end
res = table(solveIds, strings(nS,1), zeros(nS,1), zeros(nS,1), zeros(nS,1), zeros(nS,1), ...
            'VariableNames', {'id', 'location', 'Vcontact_V', 'R_kOhm', 'spread_pct', 'conservErr_pct'});
fprintf('\nSolving %d contacts (unit current %g A, idle contacts %s):\n', nS, I0, p.electrode.idleContacts);
for j = 1:nS
    s = solveIds(j);   row = find(sites.id == s);
    es = find(elemSite == s);
    K = K0;
    if p.electrode.idleContacts == "insulating"           % this contact becomes the conductor
        K = K + assembleConductionP1(tl.nodes_m, tl.tets(:, es), ...
                                     (p.electrode.contactSigma_Spm - sIdle) * [1; 1; 1]);
    end
    Te = tl.tets(:, es);
    F = accumarray(Te(:), reshape(repmat(I0 / sum(vol(es)) * vol(es) / 4, 4, 1), [], 1), [nAll 1]);
    u = zeros(nAll, 1);   tic;
    u(free) = solveSPD(K(free, free), F(free), p.solve3d);
    Ignd = -sum(K(gnd, :) * u - F(gnd));                   % current leaving through the ground
    cn = unique(Te(:));   Vc = mean(u(cn));
    res.location(j) = where(row);   res.Vcontact_V(j) = Vc;   res.R_kOhm(j) = Vc / I0 / 1e3;
    res.spread_pct(j) = 100 * (max(u(cn)) - min(u(cn))) / Vc;
    res.conservErr_pct(j) = 100 * abs(Ignd - I0) / I0;
    V(:, j) = u;
    fprintf('   AS %2d (%s, z = %+.1f mm) %-28s R = %7.2f kOhm, contact spread %.3f%%, current check %.1e%%  [%.0f s]\n', ...
            s, string(sites.type(row)), sites.z_m(row)/mm, where(row), res.R_kOhm(j), ...
            res.spread_pct(j), res.conservErr_pct(j), toc);
end


%% ===== 5. RESULTS ====================================================================
tag = patient.id + "_3d";
plotSectionLayout(patient, geom, sites, p, fullfile(outDir, tag + "_section.png"));
tr = triangulation(tl.tets', tl.nodes_m');
for sec = unique(sites.section(ismember(sites.id, solveIds)))'
    cols = find(ismember(solveIds, sites.id(sites.section == sec)));
    plotContactResults(tr, V(:, cols), sites(ismember(sites.id, solveIds(cols)), :), patient, geom, p, ...
                       res(cols, :), fullfile(outDir, sprintf('%s_section%+d.png', tag, sec)));
end

out = struct('params', p, 'patient', patient, 'sites', sites, 'solvedIds', solveIds, 'results', res, ...
    'mesh', struct('nodes_m', tl.nodes_m, 'tets', uint32(tl.tets), 'elemCell', uint32(elemCell), ...
                   'nOriginalNodes', nN, 'dupOf', uint32(tl.dupOf)), ...
    'cells', struct('tissue', geom.cellTissue, 'fascicle', geom.cellFascicle, 'site', geom.cellSite, ...
                    'zone', geom.cellZone), ...
    'V_per_A', V, 'current_A', I0, 'created', string(datetime("now")));
save(fullfile(outDir, tag + ".mat"), 'out', '-v7.3');
fprintf('\nDone: %s.mat and figures in %s\n', tag, outDir);
disp(res);


%% ===== LOCAL FUNCTIONS ===============================================================
function x = solveSPD(A, b, opt)
%SOLVESPD Solve the symmetric positive definite system A x = b (direct or preconditioned CG).
switch opt.solver
    case "direct"
        x = A \ b;
    case "pcg"
        Lc = ichol(A, struct('type', 'ict', 'droptol', 1e-4, 'diagcomp', 1e-3));
        [x, flag, relres] = pcg(A, b, opt.pcgTol, 5000, Lc, Lc');
        if flag ~= 0, warning('pcg did not converge (flag %d, relres %.1e).', flag, relres); end
    otherwise
        error('Unknown solver "%s".', opt.solver);
end
end

function where = contactLocation(sites, fasc)
%CONTACTLOCATION Where each contact lies: fascicle(s) / epineurium crossed by its mid-line.
where = strings(height(sites), 1);
for k = 1:height(sites)
    if sites.type(k) == "cuff", where(k) = "cuff, nerve surface"; continue; end
    s = linspace(-0.5, 0.5, 201)' * sites.sideXY_m(k);
    pts = sites.centre_m(k, 1:2) + s * sites.w(k, 1:2);
    lab = zeros(size(pts, 1), 1);
    for i = 1:height(fasc)
        lab(hypot(pts(:,1) - fasc.x_m(i), pts(:,2) - fasc.y_m(i)) < fasc.rEndo_m(i)) = i;
    end
    parts = strings(0);
    for v = unique(lab)'
        name = "epineurium";   if v > 0, name = "F" + v; end
        parts(end+1) = sprintf('%s %.0f%%', name, 100 * mean(lab == v)); %#ok<AGROW>
    end
    where(k) = strjoin(parts, ", ");
end
end

function drawSection(patient, geom, siteIdsToHighlight, sites)
%DRAWSECTION Nerve, fascicles (with labels) and contact footprints, in mm.
th = linspace(0, 2*pi, 300);   R = patient.nerveRadius_m;   f = patient.fascicles;
plot(1e3*R*cos(th), 1e3*R*sin(th), 'k-', LineWidth=1);   hold on
for i = 1:height(f)
    plot(1e3*(f.x_m(i) + f.rEndo_m(i)*cos(th)), 1e3*(f.y_m(i) + f.rEndo_m(i)*sin(th)), '-', ...
         Color=[0.25 0.25 0.25], LineWidth=0.8);
    text(1e3*f.x_m(i), 1e3*f.y_m(i), "F" + i, HorizontalAlignment="center", FontWeight="bold", FontSize=8);
end
for k = 1:numel(geom.footprints)
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

function plotSectionLayout(patient, geom, sites, p, fileName)
%PLOTSECTIONLAYOUT Cross-section of the patient with all contact footprints + z layout of the model.
figure(Color="w", Position=[40 60 1400 560], Name="Section and z layout");
tl = tiledlayout(1, 3, TileSpacing="compact");
title(tl, sprintf('Patient %s (seed %d): %d fascicles, R = %.2f mm, %d inner fascicles', patient.id, ...
      patient.seed, height(patient.fascicles), patient.nerveRadius_m*1e3, patient.generatorInfo.nInner), ...
      FontWeight="bold");
nexttile;
drawSection(patient, geom, [], sites);
xlabel('x [mm]'); ylabel('y [mm]');
title({'Cross-section: fascicles and', 'contact footprints of all stations'}, FontSize=9);
nexttile([1 2]);
L = p.domain.nerveLength_m;   Lc = p.cuff.length_m;   R = patient.nerveRadius_m;   Rc = R + p.cuff.wallThickness_m;
Rm = p.domain.outerRadiusFactor * R;   zF = p.mesh3d.fineHalfLength_m;
fill(1e3*[-L L L -L]/2, 1e3*[-Rm -Rm Rm Rm], [0.88 0.94 1], EdgeColor="none"); hold on
fill(1e3*[-L L L -L]/2, 1e3*[-R -R R R], [0.98 0.90 0.75], EdgeColor="k");
for sgn = [-1 1]
    fill(1e3*[-Lc Lc Lc -Lc]/2, 1e3*sgn*[R R Rc Rc], [0.3 0.3 0.3], EdgeColor="none");
end
xline(1e3*[-zF zF], ':', Color=[0.4 0.4 0.4]);
S = sites(sites.type == "intraneural", :);
for k = 1:height(S)
    plot(1e3*S.z_m(k)*[1 1], 1e3*[-1 1]*0.05, 'r-', LineWidth=2);
end
C = sites(sites.type == "cuff", :);
plot(1e3*C.z_m, 1e3*R*ones(height(C),1), 'rs', MarkerFaceColor='r', MarkerSize=4);
for planeZ = geom.planes_m, xline(1e3*planeZ, '-', Color=[0.7 0.7 0.7 0.4]); end
xlim(1e3*[-L L]/2 * 1.02); ylim(1e3*Rm*[-1 1] * 1.05);
xlabel('z [mm]'); ylabel('radial position [mm]');
title({sprintf('z layout: nerve %.0f mm (beige), cuff %.0f mm (dark, insulating), muscle (blue)', ...
               L*1e3, Lc*1e3), 'stations (red), layer planes (grey), fine-mesh zone (dotted)'}, FontSize=9);
exportgraphics(gcf, fileName);
end

function plotContactResults(tr, V, S, patient, geom, p, res, fileName)
%PLOTCONTACTRESULTS For each contact of one station: V on the cross-section through the contact,
% V along z on the axis of every fascicle, and its second difference over one internode (activating-
% function preview), all for 1 uA.
nC = height(S);   R = patient.nerveRadius_m;   f = patient.fascicles;   nF = height(f);
uA = 1e-6;   L = p.domain.nerveLength_m;   dz = p.plot3d.afSpacing_m;
figure(Color="w", Position=[20 20 min(1900, 300*nC + 120) 1000], Name="Contacts");
tl = tiledlayout(3, nC, TileSpacing="loose", Padding="loose");
title(tl, sprintf('Patient %s, station %+d (z = %.1f mm): monopolar, 1 uA', patient.id, S.section(1), ...
      S.z_m(1)*1e3), FontWeight="bold");
n = 181;   [xg, yg] = meshgrid(linspace(-1.15, 1.15, n) * R);
z = linspace(-L/2, L/2, 2801)';
cmap = lines(nF);
for j = 1:nC
    % map at the contact plane
    Vg = interpolateP1(tr, double(V(:, j)) * uA, [xg(:), yg(:), S.z_m(j) * ones(n^2, 1)]);
    nexttile(tl, j);
    imagesc(1e3*xg(1,:), 1e3*yg(:,1), reshape(1e3*Vg, n, n)); set(gca, YDir="normal"); hold on
    colormap(gca, turbo); cb = colorbar; cb.Label.String = 'V [mV]';
    drawSection(patient, geom, S.id(j), S);
    xlim(1.15e3*R*[-1 1]); ylim(1.15e3*R*[-1 1]);
    title(sprintf('AS %d (%s), %.1f kOhm', S.id(j), string(S.type(j)), res.R_kOhm(j)), FontSize=9);
    subtitle(res.location(j), FontSize=7);
    % V(z) along the fascicle axes
    Vz = zeros(numel(z), nF);
    for i = 1:nF
        Vz(:, i) = interpolateP1(tr, double(V(:, j)) * uA, [f.x_m(i)*ones(numel(z),1), f.y_m(i)*ones(numel(z),1), z]);
    end
    nexttile(tl, nC + j);
    set(gca, ColorOrder=cmap); hold on
    plot(1e3*z, 1e3*Vz, LineWidth=1);
    xline(1e3*S.z_m(j), 'r:');  xlim(1e3*(S.z_m(j) + [-12 12]*1e-3));
    xlabel('z [mm]'); ylabel('V [mV]'); title('V on the fascicle axes', FontSize=9); grid on
    if j == nC, legend("F" + (1:nF), Location="eastoutside", FontSize=7); end
    % second difference over one internode (activating function preview, V/m^2). The solution is
    % for a positive (anodic) current; for a CATHODIC pulse of 1 uA the sign flips, and a positive
    % value means depolarisation (excitation).
    k = round(dz / (z(2) - z(1)));
    afCath = -(Vz(1:end-2*k, :) - 2*Vz(1+k:end-k, :) + Vz(1+2*k:end, :)) / (k*(z(2)-z(1)))^2;
    nexttile(tl, 2*nC + j);
    set(gca, ColorOrder=cmap); hold on
    plot(1e3*z(1+k:end-k), afCath, LineWidth=1);
    xline(1e3*S.z_m(j), 'r:');  xlim(1e3*(S.z_m(j) + [-12 12]*1e-3));
    xlabel('z [mm]'); ylabel('\Delta^2V/\Deltaz^2 [V/m^2], >0 = depol.');
    title(sprintf('cathodic \\Delta^2V, \\Deltaz = %.1f mm', dz*1e3), FontSize=9); grid on
end
exportgraphics(gcf, fileName);
end
