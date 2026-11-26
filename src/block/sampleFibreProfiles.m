function prof = sampleFibreProfiles(out, b)
%SAMPLEFIBREPROFILES FEM potential along straight fibres inside the fascicles, per unit current.
%   prof = sampleFibreProfiles(out, b)
%
%   Purpose
%     Input of the fibre model (MRG in NEURON) for the conduction-block study: for every fibre and
%     every contact, the extracellular potential V(z) along the fibre for a +1 A current injected by
%     that contact alone (return = grounded outer surface). Any stimulation pattern, polarity and
%     timing is then a weighted sum of these profiles (the FEM problem is linear).
%
%   How
%     - Fibres are straight lines parallel to z (the model is an extrusion, D-03) through the nodes of
%       the viewer grid (b.gridStepXY_m), one line every b.fibreStepXY_m, kept if the line is at least
%       b.fibreInset_m inside a fascicle circle.
%     - Each line is sampled every b.profileStepZ_m over |z| <= b.profileHalfSpan_m; every sample is
%       located in the mesh AFTER the thin-layer node duplication (out.mesh), so a point inside a
%       fascicle takes the endoneurium-side value (as in the viewer, export_viewer_data.m), and V is
%       interpolated linearly (P1).
%     - Lines that pass through a contact (a fibre cannot cross metal), leave the endoneurium, or
%       have samples not found in the mesh are excluded and counted.
%     - Mirror (b.useMirror): for every solved contact k, the contact k' with the same type, radius
%       and angle at z' = -z_k gets V_k'(x, y, z) = V_k(x, y, -z) (the z grid is symmetric, so this is
%       a flip of the profile). Exact for the extruded model with a mirror-symmetric layout (D-12).
%
%   Inputs
%     out  struct saved by main_pipeline_3d (results/main3d/<id>_3d.mat): mesh, cells, patient,
%          sites, solvedIds, V_per_A [V/A]
%     b    parameters from config/block_params.m
%
%   Output: struct
%     z_m            (1, nZ) sample positions along the fibres [m]
%     x_m, y_m       (nL, 1) fibre positions in the cross-section [m]
%     fascicle       (nL, 1) fascicle index of each fibre
%     contactId      (1, nC) contact of each profile set (solved contacts first, then mirrored ones)
%     mirroredFrom   (1, nC) solved contact a mirrored profile comes from (0 = solved directly)
%     V_perA         (nL, nZ, nC) single, potential [V] for +1 A on that contact
%     excluded       struct with counts: throughContact, offEndoneurium, notFound, kept
%     source         string, which solution was sampled (out.created)

arguments
    out (1,1) struct
    b (1,1) struct
end

%% Lateral lattice (the viewer grid) and the fibre lines -------------------------------------------
X = double(out.mesh.nodes_m);   T = double(out.mesh.tets);
elemCell = double(out.mesh.elemCell);
elemTissue = out.cells.tissue(elemCell);   elemTissue = elemTissue(:)';
patient = out.patient;   F = patient.fascicles;   nF = height(F);
ax = patient.nerveSemiAxes_m;   rotN = patient.nerveRotation_deg;

h = b.gridStepXY_m;
extX = hypot(ax(1) * cosd(rotN), ax(2) * sind(rotN));   extY = hypot(ax(1) * sind(rotN), ax(2) * cosd(rotN));
ix = -ceil(extX / h):ceil(extX / h);   iy = -ceil(extY / h):ceil(extY / h);
xs = ix * h;   ys = iy * h;
every = round(b.fibreStepXY_m / h);
assert(every >= 1 && abs(every * h - b.fibreStepXY_m) < 1e-12, ...
       'sampleFibreProfiles:step', 'fibreStepXY_m must be a multiple of gridStepXY_m.');
[GX, GY] = ndgrid(xs, ys);
[IX, IY] = ndgrid(ix, iy);
fasc2 = zeros(size(GX));
for i = 1:nF
    if isfield(patient, 'outlines')                 % Gmsh model (D-14): irregular outlines, shrunk by the inset
        P = polybuffer(polyshape(patient.outlines{i}), -b.fibreInset_m);
        inside = reshape(isinterior(P, GX(:), GY(:)), size(GX));
    else                                            % circular fascicles
        inside = hypot(GX - F.x_m(i), GY - F.y_m(i)) < F.rEndo_m(i) - b.fibreInset_m;
    end
    fasc2(inside) = i;
end
line2 = fasc2 > 0 & mod(IX, every) == 0 & mod(IY, every) == 0;
nL0 = nnz(line2);

%% Samples along z, located in the mesh ----------------------------------------------------------
nzh = round(b.profileHalfSpan_m / b.profileStepZ_m);
zs = (-nzh:nzh) * b.profileStepZ_m;                 % symmetric about z = 0 (needed by the mirror)
nZ = numel(zs);
mask = repmat(line2, 1, 1, nZ);                     % find(mask): x fastest, then y, then z
nerveEl = find(ismember(elemTissue, [1 3 5]));      % endoneurium, epineurium, contacts
[tl, bary] = locateRegularGrid(X, T(:, nerveEl), {xs, ys, zs}, mask);
tetId = nan(size(tl));   ok = ~isnan(tl);   tetId(ok) = nerveEl(tl(ok));
tetId = reshape(tetId, nL0, nZ);                    % line order = find(line2) in every z slice
bary = reshape(bary, nL0, nZ, 4);

tissueAt = nan(nL0, nZ);
tissueAt(~isnan(tetId)) = elemTissue(tetId(~isnan(tetId)));
throughContact = any(tissueAt == 5, 2);
notFound = any(isnan(tetId), 2) & ~throughContact;
offEndo = any(tissueAt ~= 1 & tissueAt ~= 5 & ~isnan(tissueAt), 2) & ~throughContact & ~notFound;
keep = ~(throughContact | notFound | offEndo);
nL = nnz(keep);
assert(nL > 0, 'sampleFibreProfiles:noFibre', 'No fibre line left after the exclusions.');

%% Potential of every solved contact along the kept lines -----------------------------------------
lineIdx = find(line2);   lineIdx = lineIdx(keep);
solved = double(out.solvedIds(:)');   nS = numel(solved);
tk = tetId(keep, :);   bk = reshape(bary(keep, :, :), [], 4);
nodesOfPt = T(:, tk(:))';                           % (nL*nZ) x 4
Vs = zeros(nL, nZ, nS, 'single');
for j = 1:nS
    u = double(out.V_per_A(:, j));
    Vs(:, :, j) = single(reshape(sum(bk .* u(nodesOfPt), 2), nL, nZ));
end

%% Mirror partners ---------------------------------------------------------------------------------
sites = out.sites;
contactId = solved;   mirroredFrom = zeros(1, nS);   Vall = Vs;
if b.useMirror
    tol = 1e-9;
    for j = 1:nS
        s = sites(sites.id == solved(j), :);
        cand = sites.type == s.type & abs(sites.z_m + s.z_m) < tol & abs(sites.r_m - s.r_m) < tol & ...
               abs(mod(sites.theta_deg - s.theta_deg + 180, 360) - 180) < 1e-6 & sites.id ~= s.id;
        partner = sites.id(cand);
        assert(numel(partner) <= 1, 'sampleFibreProfiles:mirror', 'Contact %d has %d mirror partners.', ...
               solved(j), numel(partner));
        if isempty(partner) || ismember(partner, solved), continue; end
        contactId(end+1) = partner;                 %#ok<AGROW>
        mirroredFrom(end+1) = solved(j);            %#ok<AGROW>
        Vall(:, :, end+1) = flip(Vs(:, :, j), 2);   %#ok<AGROW>  V(x, y, z) -> V(x, y, -z)
    end
end

%% Output ------------------------------------------------------------------------------------------
prof = struct();
prof.z_m = zs;
prof.x_m = GX(lineIdx);   prof.y_m = GY(lineIdx);   prof.fascicle = fasc2(lineIdx);
prof.contactId = contactId;   prof.mirroredFrom = mirroredFrom;
prof.V_perA = Vall;
prof.excluded = struct('throughContact', nnz(throughContact), 'offEndoneurium', nnz(offEndo), ...
                       'notFound', nnz(notFound), 'kept', nL);
prof.source = "results/main3d/" + b.patientId + "_3d.mat, created " + string(out.created);
end
