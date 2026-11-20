%% RUN_2D_SITES_QUICKLOOK  2D nerve cross-section with the active sites of one electrode section.
% Exploratory, QUALITATIVE model: the tutor's 2D cross-section (patient P00) plus the 6 active
% sites (AS) of one section of the electrode, solved with the PDE Toolbox, one AS at a time.
%
% What a 2D model can and cannot tell (docs/architecture.md, section 2):
%   - in 2D every AS is an infinitely long strip along z, so the current is "per metre of
%     nerve" [A/m] and the potential decays like ln(r) instead of 1/r: absolute values, and
%     therefore thresholds, are NOT meaningful;
%   - no current flows along the nerve and there is no activating function (it needs dV/dz);
%   - it IS useful to see the geometry with the electrode in place, to check that the AS
%     can be meshed, and to see qualitatively how each AS spreads current among the fascicles.
%
% Parameter sets (setting paramSet):
%   "PW6"   - (default) the group's selection, config/pw6_2d_params.m and
%             report/PW6_model_choices_and_rationale.pdf: tissue conductivities (Zelechowski 2020,
%             Kangasmaa 2025), 14 um perineurium, gold contacts 50 x 200 x 200 um, cuff contacts
%             resting on the nerve surface, contacts as ideal (equipotential) boundaries,
%             hPeri = 5 um, hElec = 10 um, station z = -10.3 mm. Fascicle layout: the tutor's
%             P00 section (group choice, see the config file).
%   "tutor" - the first quick look: the tutor's conductivities (Davis 2023), 40 um perineurium,
%             5 um contacts from sites.mat, platinum metal contacts.
%
% Three ways of modelling the AS (setting siteModel), all current-controlled:
%   "ideal" - (PW6) every AS is a hole whose boundary is an ideal contact: the active AS is held
%             at 1 V, the current that flows out of it is measured, and the solution is scaled
%             to the wanted current. Idle AS are insulating holes. Bipolar pairs are solved with
%             both AS held at fixed potentials (two solves, then the combination that gives
%             -I and +I), because with insulating idle AS the monopolar solutions do not add up.
%   "metal" - every AS is a face of metal (platinum or gold, setting metalName). The active AS
%             is held at 1 V on its boundary (a metal is an equipotential), the current that
%             flows out of it is measured on a circle around the nerve, and the solution is
%             scaled to the wanted current (the problem is linear). The idle AS are floating
%             metal pieces: they carry no net current but short-circuit the tissue they touch.
%             The PDE Toolbox ignores faceLoad(CurrentDensity) in 2D dcConduction (tested),
%             so the current cannot be injected directly inside the metal.
%   "void"  - every AS is a hole ("void", D-e proposal) with a UNIFORM current density on the
%             whole boundary of the active AS; idle AS are insulating holes.
% The metal conductivity is ~10^7-10^8 times that of the tissues, so platinum and gold give
% the same field: in this resistive model the material matters only through "metal =
% equipotential". Where the material really matters is the electrode-tissue interface
% (charge-injection capacity, safety check), which is not part of the volume conductor.
%
% Stimulation (setting bipolarPairs): monopolar = one AS against the grounded outer boundary;
% bipolar = cathode -I on one AS and anode +I on another of the same section. In "metal" and
% "void" the bipolar solutions are built by superposition of the monopolar unit-current ones,
% V_bip = -V_cathode + V_anode (linearity; in "metal" a floating AS carries zero net current,
% so it does not disturb the sum); in "ideal" they are solved directly (see above). Only pairs
% inside one section are possible in 2D: a pair along the nerve (two sections) needs 3D.
%
% Other choices made for this quick look (not decisions of the project, see docs/decisions.md):
%   - both faces of an AS active (tutor's question 3 still open);
%   - outer boundary of the surrounding medium grounded, V = 0 (D-d proposal);
%   - only the tissue sigmaT is used: in a cross-section the current flows across the fibres;
%   - shank not modelled (size unknown, tutor's question 2);
%   - CENTRAL FASCICLE ENLARGED (D-g): with the tutor's radius (0.36 + 0.04 mm) the inner AS,
%     which span r = 0.2-0.4 mm, end exactly on the outer perineurium circle (r = 0.40 mm),
%     a degenerate geometry. Here its endoneurium radius is 0.42 mm, so both inner AS sit
%     inside the endoneurium (as in a TIME, where sites are meant to be inside fascicles);
%   - cuff AS: "tutor" set centred ON the nerve surface (r = R), so that the flat face does
%     not touch the circle at a single point; "PW6" set resting on it (r = R + t/2), as chosen
%     by the group (see the tangency check in section 4).
%
% Run from the project root:   setup_paths; run_2d_sites_quicklook
% Settings can be overridden before the call, e.g.
%   setup_paths; paramSet = "tutor"; run_2d_sites_quicklook
%   setup_paths; sectionId = -3; siteModel = "void"; run_2d_sites_quicklook
%   setup_paths; bipolarPairs = [2 1]; run_2d_sites_quicklook
% Outputs: figures and a .mat file in results/2d_sites/ (file names start with the set name).

clearvars -except paramSet sectionId siteModel metalName bipolarPairs; close all; clc;
mm = 1e-3;  um = 1e-6;


%% ---------- 1. SETTINGS OF THIS QUICK LOOK ----------
if ~exist('paramSet', 'var'), paramSet = "PW6"; end        % "PW6" | "tutor" (see header)
paramSet = string(paramSet);
if paramSet == "PW6"
    params = pw6_2d_params();
    defaultSection   = [];                                  % from params.electrode.zCut_m
    defaultSiteModel = params.electrode.siteModel;          % "ideal"
    defaultMetal     = params.electrode.material;           % "gold"
    tPeri_set        = params.geometry.tPeri_m;             % 14 um
    centralREndo_m   = params.geometry.centralREndo_m;      % 0.42 mm (P00 layout)
    thicknessOverride_m = params.electrode.thicknessOverride_m;   % 50 um
    cuffPlacement    = params.electrode.cuffPlacement;      % "surface": centre at R + t/2
    hSite_m = params.mesh.hSite_m;   hMin_m = params.mesh.hMin_m;
else
    params = default_params();
    defaultSection   = -4;
    defaultSiteModel = "metal";
    defaultMetal     = "platinum";
    tPeri_set        = params.geometry.tutorSection.tPeri_m;   % 40 um
    centralREndo_m   = 0.42 * mm;  % endoneurium radius of the central fascicle (see header, D-g)
    thicknessOverride_m = [];      % keep the 5 um of sites.mat
    cuffPlacement    = "centred";  % centre on r = R (see header)
    hSite_m = 5 * um;   hMin_m = 2 * um;
end
allSites = loadActiveSites(params.electrode.sitesFile, params.geometry.tutorSection.nerveRadius_m);
if ~exist('sectionId', 'var')
    if isempty(defaultSection)                              % PW6: section at z = zCut
        sectionId = unique(allSites.section(abs(allSites.z_m - params.electrode.zCut_m) < 1e-9));
        assert(isscalar(sectionId), 'No single electrode section at z = %.2f mm.', 1e3*params.electrode.zCut_m);
    else
        sectionId = defaultSection;  % electrode section to draw: -4...-1 (z < 0) or +1...+4.
    end                              % +k is the mirror image of -k, so -4...-1 cover the 4
end                                  % orientations (shank at theta = 0, 45, 90, 135 deg).
if ~exist('siteModel', 'var'), siteModel = defaultSiteModel; end   % "ideal" | "metal" | "void"
if ~exist('metalName', 'var'), metalName = defaultMetal; end
siteModel = string(siteModel);   metalName = string(metalName);
% Bulk electrical conductivity of the AS metal [S/m] at 20 C: platinum rho = 1.06e-7 Ohm*m,
% gold rho = 2.2e-8 Ohm*m (CRC Handbook of Chemistry and Physics) [CHECK the edition before
% quoting]. Any value >> tissue (~0.1 S/m) gives the same result.
sigmaMetal_Spm = struct('platinum', 1/1.06e-7, 'gold', 1/2.2e-8);
if ~exist('bipolarPairs', 'var')
    % One row per bipolar configuration: [cathode anode], as the POSITION of the AS inside the
    % section (1-4 = shank sites in file order, 5-6 = cuff sites). For section -4:
    % 1 = r 0.9 mm at +x, 2 = r 0.3 mm at +x, 3 = r 0.3 mm at -x, 4 = r 0.9 mm at -x,
    % 5 = cuff at +y, 6 = cuff at -y.
    bipolarPairs = [2 3      % both AS inside the central fascicle, across the centre (0.6 mm)
                    2 1      % inner AS (in a fascicle) vs outer AS on the same side (0.6 mm)
                    5 6];    % the two cuff AS, across the whole nerve
end
Iline_Apm     = 1;           % injected current per metre of nerve [A/m] (2D "unit current";
                             % the problem is linear, so any other value just scales V and E)

outDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'results', '2d_sites');
if ~exist(outDir, 'dir'), mkdir(outDir); end


%% ---------- 2. NERVE GEOMETRY (patient P00 = tutor's section) ----------
R     = params.geometry.tutorSection.nerveRadius_m;       % nerve radius [m]
tPeri = tPeri_set;                                        % perineurium thickness [m]
fasc  = params.geometry.tutorSection.fascicles_m;         % [x y rEndo] per fascicle [m]
fasc(1,3) = centralREndo_m;                               % D-g workaround (see header)
nFasc = size(fasc, 1);
Rmed  = params.domain.outerRadiusFactor * R;              % outer radius of the medium [m]

tissueName  = params.materials.tissueName;                % 1 endo, 2 peri, 3 epi, 4 medium
sigmaT      = params.materials.sigmaT_Spm;                % [S/m]
tissueColor = [0.95 0.70 0.75; 0.45 0.25 0.35; 0.98 0.90 0.75; 0.85 0.93 1.00];


%% ---------- 3. ACTIVE SITES OF THE CHOSEN SECTION ----------
S  = allSites(allSites.section == sectionId, :);
nS = height(S);
assert(nS > 0, 'No active site in section %d.', sectionId);

if ~isempty(thicknessOverride_m)                          % PW6: 50 um instead of the 5 um of
    S.thickness_m(:) = thicknessOverride_m;               % sites.mat (positions unchanged)
end
isCuff = S.type == "cuff";
if cuffPlacement == "surface"                             % face resting on the nerve surface
    S.r_m(isCuff) = R + S.thickness_m(isCuff) / 2;
else                                                      % "centred": centre on r = R (header)
    S.r_m(isCuff) = R;
end
S.centre_m(isCuff, 1:2) = S.r_m(isCuff) .* [cosd(S.theta_deg(isCuff)), sind(S.theta_deg(isCuff))];

fprintf('Section %+d (z = %.1f mm): %d active sites\n', sectionId, 1e3*S.z_m(1), nS);
disp(S(:, {'id', 'type', 'r_m', 'theta_deg', 'phi_deg', 'thickness_m', 'sideXY_m'}));

% Corners of each AS rectangle in the x-y plane: centre +/- (t/2) n +/- (s/2) w.
% (n, w) is a right-handed pair, so this corner order is counterclockwise.
cornerSign = [-1 -1; 1 -1; 1 1; -1 1];                    % [along n, along w]
siteCorners = zeros(4, 2, nS);                            % corner, [x y], site
for k = 1:nS
    c = S.centre_m(k, 1:2);  n = S.n(k, 1:2);  w = S.w(k, 1:2);
    a = S.thickness_m(k) / 2;  b = S.sideXY_m(k) / 2;
    siteCorners(:, :, k) = c + cornerSign(:,1) * a * n + cornerSign(:,2) * b * w;
end

% Which tissue does each AS sit in? Sample its mid-line and classify each point.
% Needed to know if an AS crosses a perineurium (D-g).
fprintf('Tissue along each AS (fraction of its length):\n');
for k = 1:nS
    s  = linspace(-0.5, 0.5, 201)' * S.sideXY_m(k);
    p  = S.centre_m(k, 1:2) + s * S.w(k, 1:2);
    [tis, fid] = tissueAt(p, fasc, tPeri, R);
    txt = strings(0);
    for t = unique(tis)'
        if mean(tis == t) < 0.01, continue; end           % skip single sample points
        f = unique(fid(tis == t & fid > 0))';
        lab = tissueName(t);
        if ~isempty(f), lab = lab + " (fascicle " + strjoin(string(f), ",") + ")"; end
        txt(end+1) = sprintf('%s %.0f%%', lab, 100 * mean(tis == t)); %#ok<SAGROW>
    end
    fprintf('   AS %2d (%-11s): %s\n', S.id(k), string(S.type(k)), strjoin(txt, ', '));
end


%% ---------- 4. 2D GEOMETRY: circles (tissues) and rectangles (AS) ----------
% Same construction as the tutor's script (section 4), with 10-row columns so that circles
% (type 1) and 4-vertex polygons (type 2) fit in the same gd matrix.
% "metal": the rectangles are added (union) and become metal faces;
% "void":  the rectangles are subtracted and become holes.
nCirc = 2 + 2*nFasc;
gd    = zeros(10, nCirc + nS);
names = strings(1, nCirc + nS);
gd(1:4, 1) = [1; 0; 0; Rmed];   names(1) = "M";
gd(1:4, 2) = [1; 0; 0; R];      names(2) = "N";
for i = 1:nFasc
    gd(1:4, 1+2*i) = [1; fasc(i,1); fasc(i,2); fasc(i,3) + tPeri];   names(1+2*i) = "P" + i;
    gd(1:4, 2+2*i) = [1; fasc(i,1); fasc(i,2); fasc(i,3)];           names(2+2*i) = "E" + i;
end
for k = 1:nS
    gd(:, nCirc + k) = [2; 4; siteCorners(:,1,k); siteCorners(:,2,k)];
    names(nCirc + k) = "S" + k;
end
ns = char(names(:))';
if siteModel == "metal"
    sf = char(strjoin(names, "+"));
else
    sf = char("(" + strjoin(names(1:nCirc), "+") + ")-(" + strjoin(names(nCirc+1:end), "+") + ")");
end
[dl, bt] = decsg(gd, sf, ns);

g = fegeometry(dl);
nFaces = g.NumFaces;
assert(size(bt, 1) == nFaces, 'Unexpected number of faces.');
assert(g.NumEdges == size(dl, 2), 'Edge numbering differs from the dl columns.');
fprintf('\nGeometry: %d faces, %d edges\n', nFaces, g.NumEdges);

% Tissue of each face (tutor's section 6). Code 5 = metal of an AS ("metal" model only;
% with "void" no face lies inside a rectangle). An AS crossing a tissue boundary is split
% into several metal faces, all with the same faceSite.
colPeri = 3:2:nCirc;   colEndo = 4:2:nCirc;   colSite = nCirc + (1:nS);
faceTissue = zeros(nFaces, 1);   faceFasc = zeros(nFaces, 1);   faceSite = zeros(nFaces, 1);
for f = 1:nFaces
    if any(bt(f, colSite)),     faceTissue(f) = 5;  faceSite(f) = find(bt(f, colSite));
    elseif any(bt(f, colEndo)), faceTissue(f) = 1;  faceFasc(f) = find(bt(f, colEndo));
    elseif any(bt(f, colPeri)), faceTissue(f) = 2;  faceFasc(f) = find(bt(f, colPeri));
    elseif bt(f, 2),            faceTissue(f) = 3;
    else,                       faceTissue(f) = 4;
    end
end

% Edges of each AS: straight segments of dl whose midpoint lies on the AS rectangle
% (boundary of the hole, or interface metal/tissue). An AS that crosses a tissue boundary
% has some sides split into several edges.
isLine = dl(1,:) == 2;
mid    = [(dl(2,:) + dl(3,:)) / 2; (dl(4,:) + dl(5,:)) / 2];
segLen = hypot(dl(3,:) - dl(2,:), dl(5,:) - dl(4,:));
siteEdges = cell(nS, 1);   sitePerim_m = zeros(nS, 1);
for k = 1:nS
    d = mid - S.centre_m(k, 1:2)';
    a = S.n(k, 1:2) * d;   b = S.w(k, 1:2) * d;           % local coordinates of each midpoint
    ha = S.thickness_m(k) / 2;   hb = S.sideXY_m(k) / 2;   tol = 1e-3 * ha;
    onSide = (abs(abs(a) - ha) < tol & abs(b) <= hb + tol) | ...
             (abs(abs(b) - hb) < tol & abs(a) <= ha + tol);
    siteEdges{k}  = find(isLine & onSide);
    sitePerim_m(k) = sum(segLen(siteEdges{k}));
    assert(abs(sitePerim_m(k) - 2*(2*ha + 2*hb)) < 1e-9, ...
           'AS %d: boundary length %.6g m, expected %.6g m.', S.id(k), sitePerim_m(k), 2*(2*ha+2*hb));
end

% Outer boundary: arcs of the medium circle
outerEdges = find(dl(1,:) == 1 & abs(hypot(dl(2,:), dl(4,:)) - Rmed) < 1e-9);


%% ---------- 5. MODEL, MATERIALS, BOUNDARY CONDITIONS, MESH ----------
model = femodel(AnalysisType = "dcConduction", Geometry = g);
for t = 1:4
    model.MaterialProperties(find(faceTissue == t)') = ...
        materialProperties(ElectricalConductivity = sigmaT(t));
end
if siteModel == "metal"
    model.MaterialProperties(find(faceTissue == 5)') = ...
        materialProperties(ElectricalConductivity = sigmaMetal_Spm.(metalName));
end
model.EdgeBC(outerEdges) = edgeBC(Voltage = 0);           % grounded outer boundary

hTissue = params.mesh.h_m;                                % endo, peri, epi, medium [m]
Hface = {find(faceTissue == 1)', hTissue(1), find(faceTissue == 2)', hTissue(2), ...
         find(faceTissue == 3)', hTissue(3)};
if siteModel == "metal", Hface = [Hface, {find(faceTissue == 5)', hSite_m}]; end
model = generateMesh(model, Hmax = hTissue(4), Hmin = hMin_m, Hgrad = params.mesh.hGrad, ...
                     Hface = Hface, Hedge = {[siteEdges{:}], hSite_m});
msh = model.Mesh;
fprintf('Mesh: %d nodes, %d triangles (quadratic)\n', size(msh.Nodes, 2), size(msh.Elements, 2));

% Check that the AS edges found from dl are the right geometry edges: their mesh nodes
% must lie on the AS rectangle.
for k = 1:nS
    nd = findNodes(msh, "region", "Edge", siteEdges{k});
    d  = msh.Nodes(:, nd) - S.centre_m(k, 1:2)';
    a  = abs(S.n(k, 1:2) * d);   b = abs(S.w(k, 1:2) * d);
    ok = all(a <= S.thickness_m(k)/2 * (1 + 1e-6) & b <= S.sideXY_m(k)/2 * (1 + 1e-6));
    assert(ok, 'AS %d: its edges are not on the AS rectangle (edge numbering problem).', S.id(k));
end


%% ---------- 6. SOLVE: monopolar, one problem per AS ----------
% Every monopolar solution is scaled to a current Iline_Apm LEAVING the AS into the tissue
% (anodic; a cathodic pulse is the same solution with -V). The current is measured as the
% flux of J through a circle of radius 1.3 R, which encloses all the AS: idle metal AS carry
% no net current, so the flux equals the current of the active AS.
%   "metal", "ideal": active AS boundary at 1 V -> measure the current -> scale.
%   "void" : uniform current density on the active AS boundary = current / perimeter
%            [A/m^2]; the measured flux is only a check (it must equal Iline_Apm).
% In "ideal" the current of the active AS is also computed as the "reaction" of the FE system
% on its nodes, sum(K(nodes,:) * V) with K the stiffness matrix (the method used for the
% bipolar pairs; checked on a hole in a disk: 0.003% from the analytical value). The two
% values must agree.
isFixedV = siteModel ~= "void";                           % active AS held at a fixed potential
siteNodes = cell(nS, 1);
for k = 1:nS, siteNodes{k} = findNodes(msh, "region", "Edge", siteEdges{k}); end
if siteModel == "ideal"
    FE = assembleFEMatrices(model, "none");               % stiffness matrix without BCs
    siteCurrent = @(u, k) full(sum(FE.K(siteNodes{k}, :) * u));   % [A/m] leaving AS k
end
nN = size(msh.Nodes, 2);
V  = zeros(nN, nS);   Ex = zeros(nN, nS);   Ey = zeros(nN, nS);   % [V], [V/m]
Iraw_Apm = zeros(nS, 1);   Ibox_Apm = nan(nS, 1);         % flux before scaling [A/m]
rC  = 1.3 * R;   thC = linspace(0, 2*pi, 1441);   thC(end) = [];
for k = 1:nS
    mk = model;                                           % femodel is a value object: a copy
    if isFixedV
        mk.EdgeBC(siteEdges{k}) = edgeBC(Voltage = 1);
    else
        mk.EdgeLoad(siteEdges{k}) = edgeLoad(SurfaceCurrentDensity = Iline_Apm / sitePerim_m(k));
    end
    res = solve(mk);
    Iraw_Apm(k) = circleFlux(res, rC, thC);
    if siteModel == "ideal", Ibox_Apm(k) = siteCurrent(res.ElectricPotential, k); end
    if isFixedV, scale = Iline_Apm / Iraw_Apm(k); else, scale = 1; end
    V(:, k)  = scale * res.ElectricPotential;
    Ex(:, k) = scale * res.ElectricField.Ex;
    Ey(:, k) = scale * res.ElectricField.Ey;
end

Vsite = arrayfun(@(k) mean(V(siteNodes{k}, k)), (1:nS)');  % potential of the active AS [V]

fprintf('\nMonopolar, %g A/m leaving the AS (site model "%s"%s):\n', Iline_Apm, siteModel, ...
        string(repmat(", " + metalName, siteModel == "metal")));
for k = 1:nS
    if siteModel == "ideal"
        fprintf('   AS %2d (%-11s): V_AS = %8.3f V   (1 V gave %.4f A/m; reaction %.4f A/m)\n', ...
                S.id(k), string(S.type(k)), Vsite(k), Iraw_Apm(k), Ibox_Apm(k));
    elseif siteModel == "metal"
        fprintf('   AS %2d (%-11s): V_AS = %8.3f V   (1 V gave %.4f A/m)\n', ...
                S.id(k), string(S.type(k)), Vsite(k), Iraw_Apm(k));
    else
        fprintf('   AS %2d (%-11s): V_AS = %8.3f V   (flux check: %.4f A/m)\n', ...
                S.id(k), string(S.type(k)), Vsite(k), Iraw_Apm(k));
    end
end
assert(all(Vsite > 0), 'Unexpected sign of the AS potential: check the current convention.');
if siteModel == "ideal"
    fprintf('Reaction vs circle flux: max relative difference %.1e\n', ...
            max(abs(Ibox_Apm - Iraw_Apm) ./ Iraw_Apm));
end

% "metal" check: an idle AS must be (nearly) an equipotential: spread of V on its boundary
% relative to the largest potential difference in the nerve.
if siteModel == "metal"
    spread = 0;
    for k = 1:nS
        for j = setdiff(1:nS, k)
            vj = V(siteNodes{j}, k);                      % (range() needs a toolbox)
            spread = max(spread, (max(vj) - min(vj)) / (max(V(:,k)) - min(V(:,k))));
        end
    end
    fprintf('Idle metal AS: max potential spread on an AS = %.1e of the range of V\n', spread);
end

% Linearity check: doubling the excitation must double V
mk = model;
if isFixedV, mk.EdgeBC(siteEdges{1}) = edgeBC(Voltage = 2);
else, mk.EdgeLoad(siteEdges{1}) = edgeLoad(SurfaceCurrentDensity = 2*Iline_Apm/sitePerim_m(1)); end
res2 = solve(mk);
V1 = V(:,1);                                              % unscaled solution of AS 1
if isFixedV, V1 = V1 * Iraw_Apm(1) / Iline_Apm; end
fprintf('Linearity check (2x excitation): max |V2 - 2 V1| / max|V1| = %.1e\n', ...
        max(abs(res2.ElectricPotential - 2*V1)) / max(abs(V1)));


%% ---------- 6b. BIPOLAR: cathode c (current -I, drawn out of the tissue), anode a (+I) ----------
% "metal", "void": superposition, V = -V_c + V_a.
% "ideal": two solves with both AS at fixed potentials, (V_c, V_a) = (1, 0) and (0, 1); the
% current of each AS is the FE reaction on its nodes, and the combination
% alpha*sol1 + beta*sol2 is chosen so that I_c = -I and I_a = +I. Both AS stay equipotential.
nB = size(bipolarPairs, 1);
assert(all(bipolarPairs(:) >= 1 & bipolarPairs(:) <= nS), 'bipolarPairs: positions 1...%d.', nS);
Vb = zeros(nN, nB);   Exb = zeros(nN, nB);   Eyb = zeros(nN, nB);
bipLabel = strings(nB, 1);
fprintf('\nBipolar, %g A/m (cathode -, anode +):\n', Iline_Apm);
for b = 1:nB
    c = bipolarPairs(b, 1);   a = bipolarPairs(b, 2);
    if siteModel == "ideal"
        sol = cell(2, 1);   Icur = zeros(2, 2);           % Icur(site, solve): rows c, a
        for q = 1:2
            mk = model;
            mk.EdgeBC(siteEdges{c}) = edgeBC(Voltage = double(q == 1));
            mk.EdgeBC(siteEdges{a}) = edgeBC(Voltage = double(q == 2));
            sol{q} = solve(mk);
            Icur(:, q) = [siteCurrent(sol{q}.ElectricPotential, c); ...
                          siteCurrent(sol{q}.ElectricPotential, a)];
        end
        coef = Icur \ [-Iline_Apm; Iline_Apm];            % [alpha; beta]
        Vb(:, b)  = coef(1) * sol{1}.ElectricPotential  + coef(2) * sol{2}.ElectricPotential;
        Exb(:, b) = coef(1) * sol{1}.ElectricField.Ex   + coef(2) * sol{2}.ElectricField.Ex;
        Eyb(:, b) = coef(1) * sol{1}.ElectricField.Ey   + coef(2) * sol{2}.ElectricField.Ey;
        Inet = coef(1) * circleFlux(sol{1}, rC, thC) + coef(2) * circleFlux(sol{2}, rC, thC);
        fprintf('   (net current to the outer boundary: %.1e A/m, should be ~0)\n', Inet);
    else
        Vb(:, b)  = -V(:, c)  + V(:, a);
        Exb(:, b) = -Ex(:, c) + Ex(:, a);
        Eyb(:, b) = -Ey(:, c) + Ey(:, a);
    end
    bipLabel(b) = sprintf('AS %d(-) / AS %d(+)', S.id(c), S.id(a));
    fprintf('   %-20s V_cathode = %8.3f V, V_anode = %8.3f V, voltage between them = %.3f V\n', ...
            bipLabel(b), mean(Vb(siteNodes{c}, b)), mean(Vb(siteNodes{a}, b)), ...
            mean(Vb(siteNodes{a}, b)) - mean(Vb(siteNodes{c}, b)));
end


%% ---------- 6c. FIELD IN EACH FASCICLE, for every configuration ----------
endoNodes = cell(nFasc, 1);
for i = 1:nFasc
    endoNodes{i} = findNodes(msh, "region", "Face", find(faceTissue == 1 & faceFasc == i)');
end
siteLabel = "AS " + S.id + " (" + string(S.type) + ")";
cfgLabel  = [siteLabel; bipLabel];                        % monopolar first, then bipolar
Emag      = hypot([Ex, Exb], [Ey, Eyb]);                  % |E| at the nodes [V/m]
meanEendo = zeros(nS + nB, nFasc);
for q = 1:nS + nB
    for i = 1:nFasc, meanEendo(q, i) = mean(Emag(endoNodes{i}, q)); end
end
fascLabel = "F" + (1:nFasc);
fprintf('\nMean |E| in the endoneurium of each fascicle [V/m per A/m]:\n');
disp(array2table(round(meanEendo, 1), 'RowNames', cellstr(cfgLabel), ...
                 'VariableNames', cellstr(fascLabel)));


%% ---------- 7. PLOTS ----------
th = linspace(0, 2*pi, 200);
xyNodes_mm = 1e3 * msh.Nodes;
tri = msh.Elements(1:3, :)';                              % corner nodes only, for patch plots

% 7a. Cross-section with the AS
figure(Color="w", Position=[40 80 1300 520], Name="Cross-section with AS");
tiledlayout(1, 2, TileSpacing="compact");
nexttile
fill(1e3*R*1.25*cos(th), 1e3*R*1.25*sin(th), tissueColor(4,:), EdgeColor="none"); hold on
fill(1e3*R*cos(th), 1e3*R*sin(th), tissueColor(3,:), EdgeColor=[0.6 0.45 0.2], LineWidth=1.5);
for i = 1:nFasc
    rO = fasc(i,3) + tPeri;
    fill(1e3*(fasc(i,1) + rO*cos(th)), 1e3*(fasc(i,2) + rO*sin(th)), tissueColor(2,:), EdgeColor="none");
    fill(1e3*(fasc(i,1) + fasc(i,3)*cos(th)), 1e3*(fasc(i,2) + fasc(i,3)*sin(th)), tissueColor(1,:), EdgeColor="none");
    text(1e3*fasc(i,1), 1e3*fasc(i,2), "F" + i, HorizontalAlignment="center", FontWeight="bold");
end
for k = 1:nS
    plot(1e3*siteCorners([1:4 1],1,k), 1e3*siteCorners([1:4 1],2,k), 'r-', LineWidth=3);
    lp = 1e3 * (S.centre_m(k,1:2) + 0.12*mm * S.n(k,1:2) .* (1 - 2*isCuff(k)));
    text(lp(1), lp(2), string(S.id(k)), Color="r", FontWeight="bold", HorizontalAlignment="center");
end
axis equal; box on; xlim(1.2e3*R*[-1 1]); ylim(1.2e3*R*[-1 1]);
xlabel('x [mm]'); ylabel('y [mm]');
title(sprintf('Section %+d: AS (red) on the cross-section', sectionId));

nexttile                                                  % mesh around the first inner AS
kZoom = find(~isCuff & S.r_m < 0.5*mm, 1);  if isempty(kZoom), kZoom = 1; end
patch(Faces=tri, Vertices=xyNodes_mm', FaceColor="w", EdgeColor=[0.3 0.3 0.3], LineWidth=0.2);
hold on; drawOutlines(R, fasc, tPeri, siteCorners);
cz = 1e3 * S.centre_m(kZoom, 1:2);
axis equal; box on; xlim(cz(1) + 0.18*[-1 1]); ylim(cz(2) + 0.18*[-1 1]);
xlabel('x [mm]'); ylabel('y [mm]');
title(sprintf('Mesh around AS %d', S.id(kZoom)));
modelTag = siteModel;  if siteModel == "metal", modelTag = metalName; end
if siteModel == "ideal", modelTag = "ideal " + metalName; end
fileTag  = sprintf('%s_section%+d_%s', paramSet, sectionId, strrep(modelTag, " ", "-"));  % e.g. PW6_section-4_ideal-gold
exportgraphics(gcf, fullfile(outDir, fileTag + "_geometry.png"));

geo = struct('msh', msh, 'tri', tri, 'xy_mm', xyNodes_mm, 'R', R, 'fasc', fasc, ...
             'tPeri', tPeri, 'siteCorners', siteCorners);

% 7b. Monopolar: potential (top row = V, the two figures split the 6 AS) and |E|
plotMaps(geo, V, Emag(:, 1:nS), siteLabel, sprintf('Monopolar (%s), %g A/m', modelTag, Iline_Apm), ...
         fullfile(outDir, fileTag + "_monopolar.png"));

% 7c. Bipolar: top row V, bottom row |E|
plotMaps(geo, Vb, Emag(:, nS+1:end), bipLabel, sprintf('Bipolar (%s), %g A/m', modelTag, Iline_Apm), ...
         fullfile(outDir, fileTag + "_bipolar.png"));

% 7d. Which fascicle does each configuration reach? Mean |E| in each endoneurium,
%     normalised by the largest value of that configuration
figure(Color="w", Position=[100 100 900 520], Name="Field per fascicle");
heatmap(cellstr(fascLabel), cellstr(cfgLabel), round(meanEendo ./ max(meanEendo, [], 2), 2), ...
        Colormap=parula, ColorLimits=[0 1]);
title(sprintf('Section %+d (%s): mean |E| in each endoneurium / max over fascicles', ...
              sectionId, modelTag));
exportgraphics(gcf, fullfile(outDir, fileTag + "_fascicles.png"));

% Light theme (as in the tutor's script, section 15)
for fig = findall(0, 'Type', 'figure')'
    if exist('theme', 'file'), theme(fig, "light"); end
end

save(fullfile(outDir, fileTag + ".mat"), 'paramSet', 'params', 'sectionId', 'siteModel', 'metalName', 'S', 'fasc', ...
     'tPeri', 'R', 'Rmed', 'Iline_Apm', 'msh', 'V', 'Ex', 'Ey', 'Vsite', 'bipolarPairs', 'Vb', ...
     'Exb', 'Eyb', 'cfgLabel', 'meanEendo', 'faceTissue', 'faceFasc', 'faceSite');
fprintf('\nDone. Figures and data in %s (files %s_*)\n', outDir, fileTag);


%% ---------- local functions ----------
function I_Apm = circleFlux(res, rC, thC)
%CIRCLEFLUX Current per metre [A/m] flowing out through the circle of radius rC [m] centred
% at the origin: integral of J . n along the circle, sampled at the angles thC.
J = interpolateCurrentDensity(res, rC*cos(thC), rC*sin(thC));
Jn = J.Jx(:)' .* cos(thC) + J.Jy(:)' .* sin(thC);
I_Apm = sum(Jn) * rC * (2*pi / numel(thC));
end

function [tis, fid] = tissueAt(p, fasc, tPeri, R)
%TISSUEAT Tissue code (1 endo, 2 peri, 3 epi, 4 medium) and fascicle index of points p [m] (N x 2).
tis = 4 * ones(size(p, 1), 1);   fid = zeros(size(p, 1), 1);
tis(hypot(p(:,1), p(:,2)) <= R) = 3;
for i = 1:size(fasc, 1)
    d = hypot(p(:,1) - fasc(i,1), p(:,2) - fasc(i,2));
    inP = d <= fasc(i,3) + tPeri;   inE = d <= fasc(i,3);
    tis(inP & ~inE) = 2;   tis(inE) = 1;   fid(inP) = i;
end
end

function plotMaps(geo, Vdata, Edata, labels, figTitle, fileName)
%PLOTMAPS One column per configuration: V (top row) and log10 |E| (bottom row), nerve zoom.
nC = size(Vdata, 2);
figure(Color="w", Position=[60 40 min(1800, 330*nC + 100) 760], Name=figTitle);
tl = tiledlayout(2, nC, TileSpacing="compact");
title(tl, figTitle, FontWeight="bold");
inNerve = hypot(geo.msh.Nodes(1,:), geo.msh.Nodes(2,:)) <= 1.1*geo.R;
for q = 1:nC
    for row = 1:2
        nexttile(tl, (row-1)*nC + q);
        if row == 1
            data = Vdata(:, q);   lim = prctile(data(inNerve), [0.5 99.5]);
            cmap = turbo;  cLab = 'V [V]';  tTxt = labels(q) + ": V";
        else
            data = log10(Edata(:, q));   lim = prctile(data(inNerve), [1 99.5]);
            cmap = parula; cLab = 'log_{10}|E| [V/m]';  tTxt = labels(q) + ": |E|";
        end
        patch(Faces=geo.tri, Vertices=geo.xy_mm', FaceVertexCData=data, FaceColor="interp", ...
              EdgeColor="none");
        hold on; drawOutlines(geo.R, geo.fasc, geo.tPeri, geo.siteCorners);
        clim(lim);  colormap(gca, cmap);  cb = colorbar;  cb.Label.String = cLab;
        axis equal; box on; xlim(1.1e3*geo.R*[-1 1]); ylim(1.1e3*geo.R*[-1 1]);
        title(tTxt, FontSize=9);
    end
end
exportgraphics(gcf, fileName);
end

function drawOutlines(R, fasc, tPeri, siteCorners)
%DRAWOUTLINES Nerve and perineurium outlines (black) and AS rectangles (red), in mm.
th = linspace(0, 2*pi, 200);
plot(1e3*R*cos(th), 1e3*R*sin(th), 'k-', LineWidth=0.8);
for i = 1:size(fasc, 1)
    rO = fasc(i,3) + tPeri;
    plot(1e3*(fasc(i,1) + rO*cos(th)), 1e3*(fasc(i,2) + rO*sin(th)), 'k-', LineWidth=0.5);
end
for k = 1:size(siteCorners, 3)
    plot(1e3*siteCorners([1:4 1],1,k), 1e3*siteCorners([1:4 1],2,k), 'r-', LineWidth=2.5);
end
end
