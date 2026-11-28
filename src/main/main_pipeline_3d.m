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
% Model in short (choices and reasons: docs/decisions.md D-02 ... D-14)
%   - mesher (p.mesh3d.mesher, D-14): "gmsh" (default) = Marco's model: elliptical nerve (D-10),
%     irregular fascicles from the chosen patient (p.patient: our generator, a virtual patient or a
%     JSON file), one random lead insertion per station (the PDMS lead bends around or pierces the
%     fascicles, pierce probability by fascicle size and angle, src/insertion/), pads placed along
%     the lead, displaced fascicles, oval insulating PDMS sleeve, flat cuff pads on the oval, 5 um
%     elements on the pad edges. The mesh is saved in results/mesh/ and reused (p.mesh3d.cacheMode).
%     "pdetoolbox" = the previous model: straight shanks, circular fascicles, void cuff wall
%     (D-03, D-05, D-11);
%   - tissues: endoneurium (anisotropic, sigma_z = longitudinal), epineurium, muscle (anisotropic);
%     PERINEURIUM = thin layer on the fascicle boundary (potential jump, sheet resistance rho*t);
%   - active contact: high-conductivity cell with a uniform source I/volume (equipotential, carries
%     exactly I); idle contacts insulating (as in 2D);
%   - outer (lateral) muscle surface grounded (p.domain.lateralBC); end planes z = +/- L/2 insulated
%     by default since 2026-11-27 (p.domain.endBC, D-17: grounded ends made the cut fibre ends fire
%     under an anode in the block study; "ground" is still available).
%   The FE system is assembled by the project's own P1 code (src/fem/), checked by
%   tests/testThinLayerFEM.m, and solved by pcg (D-14) for each contact.
%
% Units: SI inside the code (m, S/m, A, V); mm, um, uA, mV only in prints and plots.
% Run: press Run in the editor, or type main_pipeline_3d (the script adds the project paths itself).
% Time: anatomy + insertions ~ 1 min; Gmsh mesh ~ minutes (only the first time: then reused);
%   each contact ~ 1 pcg solve. Needs Python with gmsh (pip install gmsh) for the Gmsh mesher.

clear; close all; clc;
mm = 1e-3;   um = 1e-6;
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));   % project root
addpath(root);  setup_paths;
outDir = fullfile(root, 'results', 'main3d');
if ~exist(outDir, 'dir'), mkdir(outDir); end


%% ===== 1-4a. PATIENT, CONTACTS, GEOMETRY, MESH, FE SYSTEM ===============================
% src/main/build3DModel.m: patient and lead insertions (gmshAnatomy), contacts (personalDevice3D),
% mesh (Gmsh, or loaded from results/mesh/), perineurium thin layer, assembly of the stiffness
% matrix. As a function it is also called by the mesh-convergence study (scripts/run_convergence_3d.m)
% and can be called in a loop over patients / toy experiments.
p = pw6_3d_params();
M = build3DModel(p);
patient = M.patient;   sites = M.sites;   geom = M.geom;   tl = M.tl;   cuffFit = M.cuffFit;
elemCell = M.elemCell;   nN = M.nN;


%% ===== 4b. SOLUTION, one contact at a time =============================================
% src/main/solve3DContacts.m: unit current on each contact. Every contact is saved as soon as it is
% solved in results/main3d/<id>_contacts/<mesh file name>/ and reloaded by the next run if the mesh
% and the parameters are the same (resume after a stop); one folder per mesh, so switching between
% meshes (patients, toy experiments) keeps the solved contacts of each. The contacts of a previous
% complete run (results/main3d/<id>_3d.mat) seed that cache if they were solved on the same mesh.
if isempty(p.solve3d.sections), solveIds = sites.id;
else, solveIds = sites.id(ismember(sites.section, p.solve3d.sections));
end
I0 = p.stim.unitCurrent_A;
[~, meshTag] = fileparts(M.meshFile);
[V, res] = solve3DContacts(M, solveIds, struct('cacheDir', fullfile(outDir, patient.id + "_contacts", meshTag), ...
                                                'prevFile', fullfile(outDir, patient.id + "_3d.mat")));


%% ===== 5. RESULTS ====================================================================
% The solution is saved BEFORE the figures, so that a slow or interrupted plot never leaves the
% downstream scripts (export_viewer_data, export_block_profiles) with an old <id>_3d.mat.
tag = patient.id + "_3d";
out = struct('params', M.p, 'patient', patient, 'sites', sites, 'cuffFit', cuffFit, ...
    'cuffAudit', M.cuffAudit, 'meshFile', string(M.meshFile), 'solvedIds', solveIds, 'results', res, ...
    'mesh', struct('nodes_m', tl.nodes_m, 'tets', uint32(tl.tets), 'elemCell', uint32(elemCell), ...
                   'nOriginalNodes', nN, 'dupOf', uint32(tl.dupOf)), ...
    'cells', struct('tissue', geom.cellTissue, 'fascicle', geom.cellFascicle, 'site', geom.cellSite, ...
                    'zone', geom.cellZone), ...
    'V_per_A', V, 'current_A', I0, 'created', string(datetime("now")));
save(fullfile(outDir, tag + ".mat"), 'out', '-v7.3');
fprintf('\nSaved %s.mat in %s\n', tag, outDir);
disp(res);

% Figures. The points where V is sampled (fascicle axes, cross-section planes) are located in the
% mesh once and reused for every contact: locating them again for each of the 48 contacts with
% pointLocation on the ~4.5 M-element mesh took hours (2026-11-26).
plotSectionLayout(patient, geom, sites, p, fullfile(outDir, tag + "_section.png"));
zAxis_m = linspace(-p.domain.nerveLength_m/2, p.domain.nerveLength_m/2, 2801)';
axisLoc = locateOnLines(tl, patient.fascicles.x_m, patient.fascicles.y_m, zAxis_m);
for sec = unique(sites.section(ismember(sites.id, solveIds)))'
    cols = find(ismember(solveIds, sites.id(sites.section == sec)));
    plotContactResults(tl, V(:, cols), sites(ismember(sites.id, solveIds(cols)), :), patient, geom, p, ...
                       res(cols, :), fullfile(outDir, sprintf('%s_section%+d.png', tag, sec)), axisLoc, zAxis_m);
    close(gcf);
    fprintf('  figure of station %+d written\n', sec);
end
fprintf('\nDone: figures in %s\n', outDir);


%% ===== LOCAL FUNCTIONS ===============================================================
function drawSection(patient, geom, siteIdsToHighlight, sites)
%DRAWSECTION Nerve, fascicles (with labels) and contact footprints, in mm.
th = linspace(0, 2*pi, 300);   f = patient.fascicles;
nb = nerveOutline(patient.nerveSemiAxes_m, patient.nerveRotation_deg, 300);
plot(1e3*nb([1:end 1], 1), 1e3*nb([1:end 1], 2), 'k-', LineWidth=1);   hold on
for i = 1:height(f)
    if isfield(patient, 'outlines')                          % Gmsh model: the real (irregular) outlines
        P = patient.outlines{i};   plot(1e3*P([1:end 1], 1), 1e3*P([1:end 1], 2), '-', Color=[0.25 0.25 0.25], LineWidth=0.8);
    else
        plot(1e3*(f.x_m(i) + f.rEndo_m(i)*cos(th)), 1e3*(f.y_m(i) + f.rEndo_m(i)*sin(th)), '-', ...
             Color=[0.25 0.25 0.25], LineWidth=0.8);
    end
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
title(tl, sprintf('Patient %s (seed %d): %d fascicles, %s %.2f x %.2f mm (semi-axes), %d inner fascicles', ...
      patient.id, patient.seed, height(patient.fascicles), patient.nerveShape, patient.nerveSemiAxes_m*1e3, ...
      patient.generatorInfo.nInner), FontWeight="bold");
nexttile;
drawSection(patient, geom, [], sites);
xlabel('x [mm]'); ylabel('y [mm]');
title({'Cross-section: fascicles and', 'contact footprints of all stations'}, FontSize=9);
nexttile([1 2]);
% side view (schematic): the nerve is drawn with its major semi-axis
L = p.domain.nerveLength_m;   Lc = p.cuff.length_m;   R = max(patient.nerveSemiAxes_m);   Rc = R + p.cuff.wallThickness_m;
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

function plotContactResults(mesh, V, S, patient, geom, p, res, fileName, axisLoc, z)
%PLOTCONTACTRESULTS For each contact of one station: V on the cross-section through the contact,
% V along z on the axis of every fascicle, and its second difference over one internode (activating-
% function preview), all for 1 uA. mesh = tl (nodes_m 3 x nNodes, tets 4 x nElem); axisLoc =
% locateOnLines on the fascicle axes at the positions z [m] (column).
nC = height(S);   R = max(patient.nerveSemiAxes_m);   f = patient.fascicles;   nF = height(f);
uA = 1e-6;   dz = p.plot3d.afSpacing_m;
figure(Color="w", Position=[20 20 min(1900, 300*nC + 120) 1000], Name="Contacts");
tl = tiledlayout(3, nC, TileSpacing="loose", Padding="loose");
title(tl, sprintf('Patient %s, station %+d (z = %.1f mm): monopolar, 1 uA', patient.id, S.section(1), ...
      S.z_m(1)*1e3), FontWeight="bold");
n = 181;   xs = linspace(-1.15, 1.15, n) * R;
cmap = lines(nF);
planeZ = [];   planeLoc = {};                    % cross-section planes located so far
for j = 1:nC
    % map at the contact plane (located once per distinct z of the station)
    iz = find(planeZ == S.z_m(j), 1);
    if isempty(iz)
        planeZ(end+1) = S.z_m(j);   planeLoc{end+1} = locateOnPlane(mesh, xs, xs, S.z_m(j));   %#ok<AGROW>
        iz = numel(planeZ);
    end
    Vg = evalP1(mesh.tets, double(V(:, j)) * uA, planeLoc{iz});
    nexttile(tl, j);
    imagesc(1e3*xs, 1e3*xs, reshape(1e3*Vg, n, n)'); set(gca, YDir="normal"); hold on   % rows = y
    colormap(gca, turbo); cb = colorbar; cb.Label.String = 'V [mV]';
    drawSection(patient, geom, S.id(j), S);
    xlim(1.15e3*R*[-1 1]); ylim(1.15e3*R*[-1 1]);
    title(sprintf('AS %d (%s), %.1f kOhm', S.id(j), string(S.type(j)), res.R_kOhm(j)), FontSize=9);
    subtitle(res.location(j), FontSize=7);
    % V(z) along the fascicle axes
    Vz = zeros(numel(z), nF);
    for i = 1:nF
        Vz(:, i) = evalP1(mesh.tets, double(V(:, j)) * uA, axisLoc(i));
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

function loc = locateOnPlane(mesh, xs_m, ys_m, z_m)
%LOCATEONPLANE Element and barycentric coordinates of the grid points (xs_m x ys_m) on the plane z = z_m,
% x fastest. locateRegularGrid on a two-plane grid (only the first plane wanted), restricted to the
% elements that cross the plane: seconds instead of minutes with pointLocation.
h = 1e-9;                                        % spacing of the dummy second plane [m]
ez = reshape(mesh.nodes_m(3, mesh.tets), 4, []);
el = find(min(ez, [], 1) <= z_m & max(ez, [], 1) >= z_m);
mask = false(numel(xs_m), numel(ys_m), 2);   mask(:, :, 1) = true;
[t, b] = locateRegularGrid(mesh.nodes_m, double(mesh.tets(:, el)), {xs_m(:)', ys_m(:)', [z_m, z_m + h]}, mask);
ok = ~isnan(t);   t(ok) = el(t(ok));
loc = struct('tet', t, 'bary', b);
end

function loc = locateOnLines(mesh, x_m, y_m, z_m)
%LOCATEONLINES Element and barycentric coordinates of the points (x_m(i), y_m(i), z_m) on lines parallel
% to the nerve axis (one struct element per line), with the same locateRegularGrid trick as
% locateOnPlane (one x, one y), restricted to the elements whose bounding box contains the line.
h = 1e-9;
ex = reshape(mesh.nodes_m(1, mesh.tets), 4, []);   ey = reshape(mesh.nodes_m(2, mesh.tets), 4, []);
xmin = min(ex, [], 1);   xmax = max(ex, [], 1);   ymin = min(ey, [], 1);   ymax = max(ey, [], 1);
clear ex ey
loc = repmat(struct('tet', [], 'bary', []), numel(x_m), 1);
mask = false(2, 2, numel(z_m));   mask(1, 1, :) = true;
for i = 1:numel(x_m)
    el = find(xmin <= x_m(i) & xmax >= x_m(i) & ymin <= y_m(i) & ymax >= y_m(i));
    [t, b] = locateRegularGrid(mesh.nodes_m, double(mesh.tets(:, el)), ...
                               {[x_m(i), x_m(i) + h], [y_m(i), y_m(i) + h], z_m(:)'}, mask);
    ok = ~isnan(t);   t(ok) = el(t(ok));
    loc(i) = struct('tet', t, 'bary', b);
end
end

function v = evalP1(tets, u, loc)
%EVALP1 P1 interpolation of the nodal field u at points already located (loc.tet, loc.bary); NaN outside.
v = nan(numel(loc.tet), 1);   ok = ~isnan(loc.tet);
v(ok) = sum(loc.bary(ok, :) .* u(tets(:, loc.tet(ok))'), 2);
end
