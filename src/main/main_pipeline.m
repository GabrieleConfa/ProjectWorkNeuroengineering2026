%% MAIN_PIPELINE  2D median-nerve model with the PW6 parameters: geometry -> mesh -> FEM solution.
%
% Purpose
%   Builds the 2D cross-section of the median nerve at one station of the electrode, meshes it
%   and solves the volume-conductor problem div(sigma grad V) = 0 with the PDE Toolbox, for
%   every contact of the station (monopolar) and for selected contact pairs (bipolar).
%
% Inputs
%   config/pw6_2d_params.m      all parameters (values and sources: report/PW6_model_choices_
%                               and_rationale.pdf; values not in PW6 are marked there)
%   data/electrode/sites.mat    contact layout provided by the tutor (via loadActiveSites)
%
% Outputs (results/main/, generated, not versioned)
%   PW6_station<k>.mat          struct "out": parameters, contacts, mesh, potentials and fields
%   PW6_station<k>_*.png        geometry and mesh, monopolar and bipolar maps, field per fascicle
%
% Model in short
%   - tissues: endoneurium, perineurium, epineurium, surrounding muscle; in 2D only the
%     conductivity ACROSS the fibres (sigmaT) is used: the current flows in the x-y plane;
%   - contacts: gold, ideal contact boundaries = holes whose boundary is an equipotential.
%     The active contact is held at 1 V, its current is read from the FE system ("reaction")
%     and the solution is scaled to the wanted current (the problem is linear). Idle contacts
%     are insulating holes;
%   - outer boundary of the muscle: grounded (V = 0);
%   - bipolar pair (cathode c, anode a): two solves with both contacts at fixed potentials,
%     then the combination that gives -I on c and +I on a, both contacts equipotential.
%
% What 2D can and cannot tell: in 2D a contact is an infinitely long strip along the nerve,
% so currents are per metre of nerve [A/m] and the potential decays like ln(r) instead of
% 1/r. Absolute values, thresholds and the activating function (second derivative of V along
% the fibres) need the 3D model; 2D gives the geometry, the mesh and a qualitative picture of
% how each contact spreads current among the fascicles.
%
% Units: SI inside the code (m, S/m, A, V); mm only in the plots.
% Run from the project root:   setup_paths; main_pipeline

clear; close all; clc;
mm = 1e-3;
root   = fileparts(fileparts(fileparts(mfilename('fullpath'))));   % project root
outDir = fullfile(root, 'results', 'main');
if ~exist(outDir, 'dir'), mkdir(outDir); end


%% ===== 1. PARAMETERS AND CONTACTS OF THE STATION =====================================
p = pw6_2d_params();

% Nerve and fascicles. Layout chosen by the group for the current runs: the tutor's P00
% section with the central fascicle enlarged (see config/pw6_2d_params.m).
assert(p.geometry.layout == "tutorP00", 'Only the "tutorP00" layout is implemented.');
R     = p.geometry.tutorSection.nerveRadius_m;             % nerve radius [m]
fasc  = p.geometry.tutorSection.fascicles_m;               % [x y rEndo] per fascicle [m]
fasc(1,3) = p.geometry.centralREndo_m;                     % central fascicle 0.42 mm
tPeri = p.geometry.tPeri_m;                                % perineurium thickness [m]
Rmed  = p.domain.outerRadiusFactor * R;                    % outer radius of the muscle [m]
nFasc = size(fasc, 1);
sigmaT = p.materials.sigmaT_Spm;                           % [S/m]: endo, peri, epi, muscle
tissueName = p.materials.tissueName;

% Contacts of the station at z = zCut
allSites = loadActiveSites(p.electrode.sitesFile, R);
station  = unique(allSites.section(abs(allSites.z_m - p.electrode.zCut_m) < 1e-9));
assert(isscalar(station), 'No single station at z = %.2f mm.', p.electrode.zCut_m / mm);
S  = allSites(allSites.section == station, :);
nS = height(S);
S.thickness_m(:) = p.electrode.thicknessOverride_m;       % 50 um (slides), not the 5 um of sites.mat
isCuff = S.type == "cuff";
assert(p.electrode.cuffPlacement == "surface", 'Only cuff placement "surface" is implemented.');
S.r_m(isCuff) = R + S.thickness_m(isCuff) / 2;             % face resting on the nerve surface
S.centre_m(isCuff, 1:2) = S.r_m(isCuff) .* [cosd(S.theta_deg(isCuff)), sind(S.theta_deg(isCuff))];
contactLabel = "AS " + S.id + " (" + string(S.type) + ")";

fprintf('Station %+d (z = %.1f mm): %d contacts, %s, %.0f um thick\n', station, ...
        S.z_m(1) / mm, nS, p.electrode.material, S.thickness_m(1) / 1e-6);
corners = contactCorners(S);                               % 4 x 2 x nS [m]
for k = 1:nS
    fprintf('   %-18s lies in: %s\n', contactLabel(k), describeTissue(S(k,:), fasc, tPeri, R, tissueName));
end


%% ===== 2. GEOMETRY ===================================================================
% Circles (muscle, nerve, perineurium and endoneurium of each fascicle) minus the contact
% rectangles (holes). decsg splits the plane into faces; bt says which shapes contain each face.
[dl, bt] = crossSectionGeometry(R, Rmed, fasc, tPeri, corners);
g = fegeometry(dl);
assert(size(bt, 1) == g.NumFaces && g.NumEdges == size(dl, 2), 'Unexpected geometry numbering.');
[faceTissue, faceFasc] = faceTissues(bt, nFasc);           % 1 endo, 2 peri, 3 epi, 4 muscle
contactEdges = findContactEdges(dl, S);                    % geometry edges of each contact
outerEdges = find(dl(1,:) == 1 & abs(hypot(dl(2,:), dl(4,:)) - Rmed) < 1e-9);
fprintf('\nGeometry: %d faces, %d edges\n', g.NumFaces, g.NumEdges);


%% ===== 3. MESH =======================================================================
model = femodel(AnalysisType = "dcConduction", Geometry = g);
for t = 1:4
    model.MaterialProperties(find(faceTissue == t)') = ...
        materialProperties(ElectricalConductivity = sigmaT(t));
end
model.EdgeBC(outerEdges) = edgeBC(Voltage = 0);            % grounded outer boundary

h = p.mesh.h_m;                                            % target size per tissue [m]
model = generateMesh(model, Hmax = h(4), Hmin = p.mesh.hMin_m, Hgrad = p.mesh.hGrad, ...
    Hface = {find(faceTissue == 1)', h(1), find(faceTissue == 2)', h(2), find(faceTissue == 3)', h(3)}, ...
    Hedge = {[contactEdges{:}], p.mesh.hSite_m});
msh = model.Mesh;
nN  = size(msh.Nodes, 2);
q   = meshQuality(msh);
fprintf('Mesh: %d nodes, %d quadratic triangles, quality mean %.2f / min %.2f\n', ...
        nN, size(msh.Elements, 2), mean(q), min(q));

% Check: the mesh nodes of each contact's edges lie on that contact's rectangle
contactNodes = cell(nS, 1);
for k = 1:nS
    contactNodes{k} = findNodes(msh, "region", "Edge", contactEdges{k});
    [a, b] = localCoords(msh.Nodes(:, contactNodes{k}), S(k,:));
    assert(all(abs(a) <= S.thickness_m(k)/2 * (1 + 1e-6) & abs(b) <= S.sideXY_m(k)/2 * (1 + 1e-6)), ...
           'Contact %d: its edges are not on its rectangle.', S.id(k));
end


%% ===== 4. FEM SOLUTION ===============================================================
I0 = p.stim.currentPerLength_Apm;                          % injected current [A/m]
FE = assembleFEMatrices(model, "none");                    % stiffness matrix, no BCs
reaction = @(u, k) full(sum(FE.K(contactNodes{k}, :) * u));   % current leaving contact k [A/m]
rFlux = 1.3 * R;                                           % check circle, encloses all contacts

% 4a. Monopolar: contact k at 1 V -> current -> scale to I0 (anodic; cathodic = -V)
mono.V  = zeros(nN, nS);   mono.Ex = zeros(nN, nS);   mono.Ey = zeros(nN, nS);
mono.Vcontact_V = zeros(nS, 1);   mono.Iat1V_Apm = zeros(nS, 1);   Icheck = zeros(nS, 1);
for k = 1:nS
    mk = model;                                            % femodel is a value object: a copy
    mk.EdgeBC(contactEdges{k}) = edgeBC(Voltage = 1);
    res = solve(mk);
    mono.Iat1V_Apm(k) = reaction(res.ElectricPotential, k);
    Icheck(k) = circleFlux(res, rFlux);
    s = I0 / mono.Iat1V_Apm(k);
    mono.V(:, k) = s * res.ElectricPotential;
    mono.Ex(:, k) = s * res.ElectricField.Ex;   mono.Ey(:, k) = s * res.ElectricField.Ey;
    mono.Vcontact_V(k) = s;                                % the contact was at 1 V
end
fprintf('\nMonopolar, %g A/m (anodic):\n', I0);
for k = 1:nS
    fprintf('   %-18s V_contact = %7.3f V   (resistance %.2f Ohm*m)\n', contactLabel(k), ...
            mono.Vcontact_V(k), 1 / mono.Iat1V_Apm(k));
end
fprintf('   check: reaction vs flux through r = 1.3 R, max difference %.2f%%\n', ...
        100 * max(abs(Icheck - mono.Iat1V_Apm) ./ mono.Iat1V_Apm));

% 4b. Bipolar pairs [cathode anode]
pairs = p.stim.bipolarPairs;   nB = size(pairs, 1);
assert(all(pairs(:) >= 1 & pairs(:) <= nS), 'bipolarPairs: positions 1...%d.', nS);
bip.V = zeros(nN, nB);   bip.Ex = zeros(nN, nB);   bip.Ey = zeros(nN, nB);
bip.Vcathode_V = zeros(nB, 1);   bip.Vanode_V = zeros(nB, 1);   bip.label = strings(nB, 1);
fprintf('\nBipolar, %g A/m (cathode -, anode +):\n', I0);
for b = 1:nB
    c = pairs(b, 1);   a = pairs(b, 2);
    sol = cell(2, 1);   Icur = zeros(2);                   % Icur(contact c/a, solve 1/2)
    for j = 1:2                                            % solve 1: c = 1 V; solve 2: a = 1 V
        mk = model;
        mk.EdgeBC(contactEdges{c}) = edgeBC(Voltage = double(j == 1));
        mk.EdgeBC(contactEdges{a}) = edgeBC(Voltage = double(j == 2));
        sol{j} = solve(mk);
        Icur(:, j) = [reaction(sol{j}.ElectricPotential, c); reaction(sol{j}.ElectricPotential, a)];
    end
    w = Icur \ [-I0; I0];                                  % weights giving -I0 on c, +I0 on a
    bip.V(:, b)  = w(1) * sol{1}.ElectricPotential + w(2) * sol{2}.ElectricPotential;
    bip.Ex(:, b) = w(1) * sol{1}.ElectricField.Ex  + w(2) * sol{2}.ElectricField.Ex;
    bip.Ey(:, b) = w(1) * sol{1}.ElectricField.Ey  + w(2) * sol{2}.ElectricField.Ey;
    bip.Vcathode_V(b) = w(1);   bip.Vanode_V(b) = w(2);    % potentials of the two contacts
    bip.label(b) = sprintf('AS %d(-) / AS %d(+)', S.id(c), S.id(a));
    Inet = w(1) * circleFlux(sol{1}, rFlux) + w(2) * circleFlux(sol{2}, rFlux);
    fprintf('   %-18s V_cathode = %7.3f V, V_anode = %7.3f V (net current to ground %.2f%%)\n', ...
            bip.label(b), w(1), w(2), 100 * Inet / I0);
end


%% ===== 5. RESULTS ====================================================================
% Mean |E| in the endoneurium of each fascicle: a qualitative measure of how much each
% configuration reaches each fascicle (no activating function in 2D).
endoNodes = arrayfun(@(i) findNodes(msh, "region", "Face", find(faceTissue == 1 & faceFasc == i)'), ...
                     1:nFasc, 'UniformOutput', false);
cfgLabel = [contactLabel; bip.label];
Emag = hypot([mono.Ex, bip.Ex], [mono.Ey, bip.Ey]);        % [V/m]
meanEendo = cell2mat(cellfun(@(nd) mean(Emag(nd, :), 1)', endoNodes, 'UniformOutput', false));
fascLabel = "F" + (1:nFasc);
fprintf('\nMean |E| in each endoneurium [V/m for %g A/m]:\n', I0);
disp(array2table(round(meanEendo), 'RowNames', cellstr(cfgLabel), 'VariableNames', cellstr(fascLabel)));

tag = sprintf('PW6_station%+d', station);
geo = struct('msh', msh, 'R', R, 'fasc', fasc, 'tPeri', tPeri, 'corners', corners);
plotGeometryAndMesh(geo, S, faceTissue, sprintf('Station %+d (z = %.1f mm)', station, S.z_m(1)/mm), ...
                    fullfile(outDir, tag + "_geometry.png"));
plotMaps(geo, mono.V, Emag(:, 1:nS), contactLabel, sprintf('Monopolar, %g A/m', I0), ...
         fullfile(outDir, tag + "_monopolar.png"));
plotMaps(geo, bip.V, Emag(:, nS+1:end), bip.label, sprintf('Bipolar, %g A/m', I0), ...
         fullfile(outDir, tag + "_bipolar.png"));
plotFascicleTable(meanEendo ./ max(meanEendo, [], 2), cfgLabel, fascLabel, ...
                  sprintf('Station %+d: mean |E| in each endoneurium / max of the row', station), ...
                  fullfile(outDir, tag + "_fascicles.png"));

out = struct('params', p, 'station', station, 'contacts', S, 'fascicles_m', fasc, ...
             'mesh', struct('nodes_m', msh.Nodes, 'elements', msh.Elements), ...
             'faceTissue', faceTissue, 'faceFascicle', faceFasc, 'current_Apm', I0, ...
             'monopolar', mono, 'bipolar', bip, 'bipolarPairs', pairs, ...
             'meanEendo_Vpm', meanEendo, 'configLabel', cfgLabel, 'created', string(datetime("now")));
save(fullfile(outDir, tag + ".mat"), 'out');
fprintf('\nDone: %s.mat and figures in %s\n', tag, outDir);


%% ===== LOCAL FUNCTIONS ===============================================================
function corners = contactCorners(S)
%CONTACTCORNERS Corners [m] of each contact rectangle in the x-y plane (4 x 2 x nContacts),
% counterclockwise: centre +/- (t/2) n +/- (s/2) w.
sgn = [-1 -1; 1 -1; 1 1; -1 1];
corners = zeros(4, 2, height(S));
for k = 1:height(S)
    corners(:, :, k) = S.centre_m(k, 1:2) + sgn(:,1) * (S.thickness_m(k)/2) * S.n(k, 1:2) ...
                                          + sgn(:,2) * (S.sideXY_m(k)/2)   * S.w(k, 1:2);
end
end

function [a, b] = localCoords(xy, site)
%LOCALCOORDS Coordinates of points xy [m] (2 x N) along the normal n and the side w of a contact.
d = xy - site.centre_m(1:2)';
a = site.n(1:2) * d;   b = site.w(1:2) * d;
end

function txt = describeTissue(site, fasc, tPeri, R, tissueName)
%DESCRIBETISSUE Tissues crossed by the mid-line of a contact, with the fraction of its length.
s = linspace(-0.5, 0.5, 201)' * site.sideXY_m;
pts = site.centre_m(1:2) + s * site.w(1:2);
tis = 4 * ones(size(pts, 1), 1);   fid = zeros(size(pts, 1), 1);
tis(hypot(pts(:,1), pts(:,2)) <= R) = 3;
for i = 1:size(fasc, 1)
    d = hypot(pts(:,1) - fasc(i,1), pts(:,2) - fasc(i,2));
    tis(d <= fasc(i,3) + tPeri & d > fasc(i,3)) = 2;   tis(d <= fasc(i,3)) = 1;
    fid(d <= fasc(i,3) + tPeri) = i;
end
parts = strings(0);
for t = unique(tis)'
    if mean(tis == t) < 0.01, continue; end
    lab = tissueName(t);
    f = unique(fid(tis == t & fid > 0))';
    if ~isempty(f), lab = lab + " (F" + strjoin(string(f), ",") + ")"; end
    parts(end+1) = sprintf('%s %.0f%%', lab, 100 * mean(tis == t)); %#ok<AGROW>
end
txt = strjoin(parts, ', ');
end

function [dl, bt] = crossSectionGeometry(R, Rmed, fasc, tPeri, corners)
%CROSSSECTIONGEOMETRY decsg geometry: circles M (muscle), N (nerve), P_i (perineurium) and
% E_i (endoneurium) united, minus the contact rectangles S_k. Columns of bt: M, N, P1, E1, ...
nF = size(fasc, 1);   nC = 2 + 2*nF;   nS = size(corners, 3);
gd = zeros(10, nC + nS);   names = strings(1, nC + nS);
gd(1:4, 1) = [1; 0; 0; Rmed];   names(1) = "M";
gd(1:4, 2) = [1; 0; 0; R];      names(2) = "N";
for i = 1:nF
    gd(1:4, 1+2*i) = [1; fasc(i,1); fasc(i,2); fasc(i,3) + tPeri];   names(1+2*i) = "P" + i;
    gd(1:4, 2+2*i) = [1; fasc(i,1); fasc(i,2); fasc(i,3)];           names(2+2*i) = "E" + i;
end
for k = 1:nS
    gd(:, nC + k) = [2; 4; corners(:,1,k); corners(:,2,k)];   names(nC + k) = "S" + k;
end
sf = "(" + strjoin(names(1:nC), "+") + ")-(" + strjoin(names(nC+1:end), "+") + ")";
[dl, bt] = decsg(gd, char(sf), char(names(:))');
end

function [faceTissue, faceFasc] = faceTissues(bt, nF)
%FACETISSUES Tissue (1 endo, 2 peri, 3 epi, 4 muscle) and fascicle index of each face, from
% the boolean table of decsg. The innermost circle is tested first.
colP = 3:2:2+2*nF;   colE = 4:2:2+2*nF;
nFaces = size(bt, 1);   faceTissue = 4 * ones(nFaces, 1);   faceFasc = zeros(nFaces, 1);
for f = 1:nFaces
    if any(bt(f, colE)),     faceTissue(f) = 1;  faceFasc(f) = find(bt(f, colE));
    elseif any(bt(f, colP)), faceTissue(f) = 2;  faceFasc(f) = find(bt(f, colP));
    elseif bt(f, 2),         faceTissue(f) = 3;
    end
end
end

function edges = findContactEdges(dl, S)
%FINDCONTACTEDGES Geometry edges of each contact: straight dl segments whose midpoint lies on
% the contact rectangle (a side crossing a tissue boundary is split into several edges).
% Checks that the edges found add up to the contact perimeter.
isLine = dl(1,:) == 2;
mid = [(dl(2,:) + dl(3,:)) / 2; (dl(4,:) + dl(5,:)) / 2];
len = hypot(dl(3,:) - dl(2,:), dl(5,:) - dl(4,:));
edges = cell(height(S), 1);
for k = 1:height(S)
    [a, b] = localCoords(mid, S(k,:));
    ha = S.thickness_m(k) / 2;   hb = S.sideXY_m(k) / 2;   tol = 1e-3 * ha;
    on = (abs(abs(a) - ha) < tol & abs(b) <= hb + tol) | (abs(abs(b) - hb) < tol & abs(a) <= ha + tol);
    edges{k} = find(isLine & on);
    assert(abs(sum(len(edges{k})) - 4*(ha + hb)) < 1e-9, 'Contact %d: wrong perimeter.', S.id(k));
end
end

function I_Apm = circleFlux(res, r)
%CIRCLEFLUX Current per metre [A/m] flowing out through the circle of radius r [m] around the
% origin (integral of J . n), used as an independent check of the contact currents.
th = linspace(0, 2*pi, 1441);   th(end) = [];
J = interpolateCurrentDensity(res, r*cos(th), r*sin(th));
I_Apm = sum(J.Jx(:)' .* cos(th) + J.Jy(:)' .* sin(th)) * r * (2*pi / numel(th));
end

function drawOutlines(geo)
%DRAWOUTLINES Nerve and perineurium outlines (black) and contacts (red), in mm.
th = linspace(0, 2*pi, 200);
plot(1e3*geo.R*cos(th), 1e3*geo.R*sin(th), 'k-', LineWidth=0.8);
for i = 1:size(geo.fasc, 1)
    rO = geo.fasc(i,3) + geo.tPeri;
    plot(1e3*(geo.fasc(i,1) + rO*cos(th)), 1e3*(geo.fasc(i,2) + rO*sin(th)), 'k-', LineWidth=0.5);
end
for k = 1:size(geo.corners, 3)
    plot(1e3*geo.corners([1:4 1],1,k), 1e3*geo.corners([1:4 1],2,k), 'r-', LineWidth=1.5);
end
end

function plotGeometryAndMesh(geo, S, faceTissue, figTitle, fileName)
%PLOTGEOMETRYANDMESH Cross-section coloured by tissue (left) and mesh around an inner contact.
col = [0.95 0.70 0.75; 0.45 0.25 0.35; 0.98 0.90 0.75; 0.85 0.93 1.00];
msh = geo.msh;   tri = msh.Elements(1:3, :)';   xy = 1e3 * msh.Nodes';
elemTissue = zeros(size(tri, 1), 1);
for f = 1:numel(faceTissue)
    elemTissue(findElements(msh, "region", "Face", f)) = faceTissue(f);
end
figure(Color="w", Position=[40 80 1300 540], Name="Geometry and mesh");
tl = tiledlayout(1, 2, TileSpacing="compact");   title(tl, figTitle, FontWeight="bold");
nexttile
patch(Faces=tri, Vertices=xy, FaceVertexCData=col(elemTissue,:), FaceColor="flat", EdgeColor="none");
hold on; drawOutlines(geo);
for k = 1:height(S)
    lp = 1e3 * (S.centre_m(k,1:2) + 0.13e-3 * S.n(k,1:2) * (1 - 2*(S.type(k) == "cuff")));
    text(lp(1), lp(2), string(S.id(k)), Color="r", FontWeight="bold", HorizontalAlignment="center");
end
for i = 1:size(geo.fasc, 1)
    text(1e3*geo.fasc(i,1), 1e3*geo.fasc(i,2) + 0.12*(i == 1), "F" + i, HorizontalAlignment="center", FontWeight="bold");
end
axis equal; box on; xlim(1.2e3*geo.R*[-1 1]); ylim(1.2e3*geo.R*[-1 1]);
xlabel('x [mm]'); ylabel('y [mm]'); title('Tissues and contacts');
nexttile
k = find(S.type ~= "cuff" & S.r_m < 0.5e-3, 1);   if isempty(k), k = 1; end
patch(Faces=tri, Vertices=xy, FaceColor="w", EdgeColor=[0.3 0.3 0.3], LineWidth=0.2);
hold on; drawOutlines(geo);
c = 1e3 * S.centre_m(k, 1:2);
axis equal; box on; xlim(c(1) + 0.2*[-1 1]); ylim(c(2) + 0.2*[-1 1]);
xlabel('x [mm]'); ylabel('y [mm]'); title(sprintf('Mesh around AS %d', S.id(k)));
exportgraphics(gcf, fileName);
end

function plotMaps(geo, Vdata, Edata, labels, figTitle, fileName)
%PLOTMAPS One column per configuration: V (top) and log10 |E| (bottom), zoom on the nerve.
msh = geo.msh;   tri = msh.Elements(1:3, :)';   xy = 1e3 * msh.Nodes';
nC = size(Vdata, 2);
inNerve = hypot(msh.Nodes(1,:), msh.Nodes(2,:)) <= 1.1 * geo.R;
figure(Color="w", Position=[60 40 min(1800, 330*nC + 100) 760], Name=figTitle);
tl = tiledlayout(2, nC, TileSpacing="compact");   title(tl, figTitle, FontWeight="bold");
for j = 1:nC
    for row = 1:2
        nexttile(tl, (row-1)*nC + j);
        if row == 1
            d = Vdata(:, j);   lim = prctile(d(inNerve), [0.5 99.5]);   cm = turbo;   cl = 'V [V]';
        else
            d = log10(Edata(:, j));   lim = prctile(d(inNerve), [1 99.5]);   cm = parula;  cl = 'log_{10}|E| [V/m]';
        end
        patch(Faces=tri, Vertices=xy, FaceVertexCData=d, FaceColor="interp", EdgeColor="none");
        hold on; drawOutlines(geo);
        clim(lim);  colormap(gca, cm);  cb = colorbar;  cb.Label.String = cl;
        axis equal; box on; xlim(1.1e3*geo.R*[-1 1]); ylim(1.1e3*geo.R*[-1 1]);
        title(labels(j) + string(repmat(": V", row == 1)) + string(repmat(": |E|", row == 2)), FontSize=9);
    end
end
exportgraphics(gcf, fileName);
end

function plotFascicleTable(M, rowLabel, colLabel, figTitle, fileName)
%PLOTFASCICLETABLE Colour table of M (rows = configurations, columns = fascicles), with values.
figure(Color="w", Position=[100 100 900 520], Name="Field per fascicle");
imagesc(M, [0 1]);  colormap(parula);  colorbar;
for r = 1:size(M, 1)
    for c = 1:size(M, 2)
        text(c, r, sprintf('%.2f', M(r, c)), HorizontalAlignment="center", ...
             Color=[1 1 1] * (M(r, c) < 0.6));
    end
end
set(gca, XTick=1:numel(colLabel), XTickLabel=colLabel, YTick=1:numel(rowLabel), ...
         YTickLabel=rowLabel, TickLength=[0 0]);
title(figTitle);
exportgraphics(gcf, fileName);
end
