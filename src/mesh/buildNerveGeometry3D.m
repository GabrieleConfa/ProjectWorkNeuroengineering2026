function geom = buildNerveGeometry3D(patient, sites, p)
%BUILDNERVEGEOMETRY3D 3D geometry of nerve + cuff + contacts + surrounding muscle.
%   geom = buildNerveGeometry3D(patient, sites, p)
%
%   How it is built (docs/decisions.md D-03, explained in docs/model_3d.md §6):
%   1. A 2D cross-section with the tissue outlines only: muscle (circle 3R), outer surface of the
%      cuff wall (circle R + wall), nerve (circle R), fascicles (circles rEndo: the perineurium is a
%      thin layer ON these circles, not a volume).
%   2. Extrusion in 5 LAYERS along z (extrude with a vector of heights): end zone, cuff zone, fine
%      zone, cuff zone, end zone. Every 2D face becomes one cell per layer.
%   3. The cuff wall ring is deleted inside the cuff (void = perfect insulator, so the nerve surface
%      under the cuff is an insulated boundary) and merged into the muscle outside the cuff.
%   4. The 48 contact boxes are inserted one by one with a 3D Boolean union that keeps the existing
%      boundaries (union(..., KeepBoundaries=[true true])): a box crossing a fascicle boundary is
%      split into one piece per tissue, and those pieces are merged into one contact cell. A cuff
%      box (pressed cuffIndent into the nerve) has a thin piece inside the nerve and a piece in the
%      deleted wall, which becomes part of the contact: the contact is embedded in the cuff wall.
%
%   Inputs
%     patient  struct from makePatient (nerveRadius_m, fascicles table)
%     sites    table from loadActiveSites, with the model thickness and the cuff-contact centres
%              already set by the caller (main_pipeline_3d.m)
%     p        parameters (config/pw6_3d_params.m)
%
%   Output: struct
%     g              fegeometry (3D, not meshed), z = -L/2 ... +L/2
%     cellTissue     (nCells,1) 1 endoneurium, 3 epineurium, 4 muscle, 5 contact
%     cellFascicle   (nCells,1) fascicle index, 0 if not endoneurium
%     cellSite       (nCells,1) site id for contact cells, 0 otherwise
%     cellZone       (nCells,1) -2/2 outside the cuff, -1/1 cuff but not fine, 0 fine zone
%     interfaceFaces faces between endoneurium and epineurium (perineurium thin layer)
%     groundFaces    faces with V = 0 (outer muscle surface; end planes if p.domain.endBC = "ground")
%     faceSize_m     (nFaces,1) target element size on each face (Inf = global maximum)
%     planes_m       z of the layer planes
%     footprints     struct array (corners_m, siteIds, kind): contact outlines in x-y, for plots
%     section2d      dl, bt, faceTissue, faceFascicle, facePoint_m
%     nMerges, time_s

arguments
    patient (1,1) struct
    sites table
    p (1,1) struct
end
t0 = tic;
R  = patient.nerveRadius_m;   F = patient.fascicles;   nF = height(F);
Rm = p.domain.outerRadiusFactor * R;
Rc = R + p.cuff.wallThickness_m;
L  = p.domain.nerveLength_m;   Lc = p.cuff.length_m;   zF = p.mesh3d.fineHalfLength_m;
delta = p.electrode.cuffIndent_m;

%% Checks on the contacts ------------------------------------------------------------------------
hz = max(sites.sideZ_m) / 2;
assert(max(abs(sites.z_m)) + hz < zF && zF < Lc/2 && Lc/2 < L/2, ...
       'Need: contacts < fine zone (%.1f mm) < cuff/2 (%.1f mm) < nerve/2 (%.1f mm).', zF*1e3, Lc/2e3, L/2e3);
cuff = sites.type == "cuff";
hw = max(sites.sideXY_m(cuff)) / 2;
assert(delta > hw^2 / (2*R), 'cuffIndent_m (%.1f um) must exceed the sagitta %.1f um.', delta*1e6, hw^2/(2*R)*1e6);
assert(all(sites.thickness_m(cuff) - delta < Rc - R), 'The cuff contacts must fit inside the cuff wall.');

%% 2D section: tissues only -------------------------------------------------------------------------
gd = zeros(10, 0);   names = strings(1, 0);
gd(:, end+1) = [1; 0; 0; Rm; zeros(6,1)];   names(end+1) = "M";
gd(:, end+1) = [1; 0; 0; Rc; zeros(6,1)];   names(end+1) = "C";
gd(:, end+1) = [1; 0; 0; R;  zeros(6,1)];   names(end+1) = "N";
for i = 1:nF
    gd(:, end+1) = [1; F.x_m(i); F.y_m(i); F.rEndo_m(i); zeros(6,1)];   names(end+1) = "E" + i; %#ok<AGROW>
end
[dl, bt] = decsg(gd, char(strjoin(names, "+")), char(names(:))');
g2 = fegeometry(dl);
nF2 = size(bt, 1);   colE = 3 + (1:nF);
faceTissue = zeros(nF2, 1);   faceFasc = zeros(nF2, 1);
for f = 1:nF2
    if any(bt(f, colE))
        i = find(bt(f, colE));   assert(isscalar(i), 'Two fascicles overlap.');
        faceTissue(f) = 1;   faceFasc(f) = i;
    elseif bt(f, 3), faceTissue(f) = 3;               % epineurium (inside the nerve)
    elseif bt(f, 2), faceTissue(f) = 6;               % cuff wall ring (R ... R + wall)
    else,            faceTissue(f) = 4;               % muscle
    end
end
facePt = interiorPoints(g2, R / 15);                  % one point safely inside every 2D face

%% Extrusion in 5 layers ----------------------------------------------------------------------------
planes = [-L/2, -Lc/2, -zF, zF, Lc/2, L/2];
H = diff(planes);   zMid = (planes(1:end-1) + planes(2:end)) / 2;   nL = numel(H);
zone = [-2 -1 0 1 2];
g = translate(extrude(g2, H), [0 0 planes(1)]);
nMerge = 0;
ring = find(faceTissue == 6);   musc = find(faceTissue == 4);
for L_ = [1 nL]                                       % outside the cuff: wall ring -> muscle
    [g, n] = mergeGroup(g, [facePt([musc; ring], :), zMid(L_) * ones(numel(musc) + numel(ring), 1)]);
    nMerge = nMerge + n;
end
voidCells = findCell(g, [facePt(ring, :), zMid(2) * ones(numel(ring),1); ...
                         facePt(ring, :), zMid(3) * ones(numel(ring),1); ...
                         facePt(ring, :), zMid(4) * ones(numel(ring),1)]);
g = deleteCell(g, unique(voidCells));                 % inside the cuff: wall ring = void

%% Contacts: one box per site, inserted with a Boolean union ----------------------------------------
nS = height(sites);
for s = 1:nS
    box = fegeometry(multicuboid(sites.thickness_m(s), sites.sideXY_m(s), sites.sideZ_m(s)));
    box = rotate(box, sites.phi_deg(s), [0 0 0], [0 0 1]);          % local x -> normal n
    box = translate(box, sites.centre_m(s, :) - [0 0 sites.sideZ_m(s)/2]);
    g = union(g, box, KeepBoundaries=[true true]);
    [g, n] = mergeGroup(g, boxPoints(sites(s, :)));                 % pieces -> one contact cell
    nMerge = nMerge + n;
end

%% Label the final cells ------------------------------------------------------------------------------
nC = g.NumCells;
cellTissue = zeros(nC, 1);   cellFascicle = zeros(nC, 1);   cellSite = zeros(nC, 1);   cellZone = nan(nC, 1);
for L_ = 1:nL
    keep = ~(faceTissue == 6 & abs(zone(L_)) < 2);    % the wall inside the cuff is void
    fl = find(keep);
    cid = findCell(g, [facePt(fl, :), zMid(L_) * ones(numel(fl), 1)]);
    assert(~any(isnan(cid)), 'A tissue reference point fell outside the geometry (layer %d).', L_);
    for q = 1:numel(fl)
        t = faceTissue(fl(q));   if t == 6, t = 4; end
        setRole(cid(q), t, faceFasc(fl(q)), 0, zone(L_));
    end
end
for s = 1:nS
    cid = unique(findCell(g, boxPoints(sites(s, :))));
    assert(isscalar(cid) && ~isnan(cid), 'Contact %d is not a single cell.', sites.id(s));
    setRole(cid, 5, 0, sites.id(s), 0);
end
assert(all(cellTissue > 0), 'Cells without a role: %s.', mat2str(find(cellTissue == 0)'));

%% Faces: thin layer, ground, mesh sizes --------------------------------------------------------------
facesOf = @(cells) unique(cell2mat(arrayfun(@(c) reshape(cellFaces(g, c), 1, []), cells(:)', ...
                                             'UniformOutput', false)));
interfaceFaces = intersect(facesOf(find(cellTissue == 1)), facesOf(find(cellTissue == 3)));
ang = deg2rad(22.5 + 45*(0:7));
latPts = [kron(ones(nL,1), Rm * [cos(ang') sin(ang')]), kron(zMid', ones(8,1))];
groundFaces = unique(nearestFace(g, latPts));
if p.domain.endBC == "ground"
    endPts = [facePt, -L/2 * ones(nF2,1); facePt, L/2 * ones(nF2,1)];
    groundFaces = unique([groundFaces(:); reshape(nearestFace(g, endPts), [], 1)]);
elseif p.domain.endBC ~= "insulated"
    error('buildNerveGeometry3D:endBC', 'Unknown endBC "%s".', p.domain.endBC);
end
faceSize = inf(g.NumFaces, 1);
nerveCells = find(ismember(cellTissue, [1 3]));
faceSize(facesOf(nerveCells(cellZone(nerveCells) ~= 0))) = p.mesh3d.hCoarse_m;
sizeRules = {find(cellTissue == 3 & cellZone == 0), p.mesh3d.hNerve_m; ...
             find(cellTissue == 1 & cellZone == 0), p.mesh3d.hEndo_m; ...
             find(cellTissue == 5),                 p.mesh3d.hContact_m};
for r = 1:size(sizeRules, 1)                          % the smallest size wins on shared faces
    fc = facesOf(sizeRules{r, 1});
    faceSize(fc) = min(faceSize(fc), sizeRules{r, 2});
end

geom = struct('g', g, 'cellTissue', cellTissue, 'cellFascicle', cellFascicle, 'cellSite', cellSite, ...
    'cellZone', cellZone, 'interfaceFaces', interfaceFaces(:), 'groundFaces', groundFaces(:), ...
    'faceSize_m', faceSize, 'planes_m', planes, 'footprints', contactFootprints(sites), ...
    'section2d', struct('dl', dl, 'bt', bt, 'faceTissue', faceTissue, 'faceFascicle', faceFasc, ...
                        'facePoint_m', facePt), ...
    'nMerges', nMerge, 'time_s', toc(t0));

    function setRole(c, t, fascIdx, siteId, z)
    %SETROLE Give cell c its role; a cell reached by two different roles is an error.
    if cellTissue(c) == 0
        cellTissue(c) = t;   cellFascicle(c) = fascIdx;   cellSite(c) = siteId;   cellZone(c) = z;
    else
        assert(cellTissue(c) == t && cellFascicle(c) == fascIdx && cellSite(c) == siteId, ...
               'Cell %d has two roles.', c);
    end
    end
end

function P = boxPoints(site)
%BOXPOINTS 3 x 5 x 3 points inside a contact box (at 80% of its half-sides), used to find all the
% pieces of the box after the union and to check its orientation.
[a, b, c] = ndgrid([-0.8 0 0.8] * site.thickness_m / 2, linspace(-0.8, 0.8, 5) * site.sideXY_m / 2, ...
                   [-0.8 0 0.8] * site.sideZ_m / 2);
P = site.centre_m + a(:) * site.n + b(:) * site.w + c(:) * [0 0 1];
end

function fp = contactFootprints(sites)
%CONTACTFOOTPRINTS Outline of every distinct contact in the x-y plane (sections +k and -k share it).
fp = struct('corners_m', {}, 'siteIds', {}, 'kind', {});
key = round([sites.centre_m(:,1:2), sites.phi_deg / 1e3] * 1e7);
[~, iu, ic] = unique(key, 'rows', 'stable');
sgn = [-1 -1; 1 -1; 1 1; -1 1];
for k = 1:numel(iu)
    s = sites(iu(k), :);
    fp(k).corners_m = s.centre_m(1:2) + sgn(:,1) * (s.thickness_m/2) * s.n(1:2) ...
                                      + sgn(:,2) * (s.sideXY_m/2) * s.w(1:2);
    fp(k).siteIds = sites.id(ic == k);
    fp(k).kind = string(s.type);
end
end

function pt = interiorPoints(g2, hmax)
%INTERIORPOINTS One point safely inside every face of a 2D geometry: the centroid of the largest
% triangle of that face in a fine 2D mesh.
m2 = generateMesh(g2, Hmax=hmax, GeometricOrder="linear");
P = m2.Mesh.Nodes;   pt = zeros(g2.NumFaces, 2);
for f = 1:g2.NumFaces
    el = findElements(m2.Mesh, "region", "Face", f);
    T = m2.Mesh.Elements(1:3, el);
    ar = abs((P(1,T(2,:)) - P(1,T(1,:))) .* (P(2,T(3,:)) - P(2,T(1,:))) - ...
             (P(1,T(3,:)) - P(1,T(1,:))) .* (P(2,T(2,:)) - P(2,T(1,:))));
    [~, b] = max(ar);
    pt(f, :) = mean(P(:, T(:, b)), 2)';
end
end

function [g, nDone] = mergeGroup(g, P)
%MERGEGROUP Merge all the cells that contain the points P (n x 3) into one, pairwise (mergeCells
% takes two adjacent cells). Cell IDs change after every merge, so they are looked up again.
nDone = 0;
while true
    ids = unique(findCell(g, P));
    ids = ids(~isnan(ids));
    if numel(ids) <= 1, return; end
    merged = false;
    for i = 1:numel(ids)-1
        for j = i+1:numel(ids)
            try
                g = mergeCells(g, [ids(i) ids(j)]);
                merged = true;   nDone = nDone + 1;
                break
            catch
            end
        end
        if merged, break; end
    end
    if ~merged
        warning('buildNerveGeometry3D:merge', 'Cells %s could not be merged (not adjacent).', mat2str(ids'));
        return
    end
end
end
