function fx = fieldExposure(out, ids, weights)
%FIELDEXPOSURE How much tissue sees a strong electric field, for 1 A of a contact pattern.
%   fx = fieldExposure(out, ids, weights)
%
%   Purpose
%     Safety check of the field strength (scripts/run_field_safety.m). A strong field electroporates
%     the cell membranes: it damages the tissue, and it raises the tissue conductivity, which our
%     linear FEM ignores. This function computes the distribution of |E| in each tissue for 1 A of
%     total pattern current. Because the problem is linear, |E| for a current I is I times it, so
%     one call answers "how much tissue exceeds field level E at current I" for every E and I
%     (exposedVolume).
%
%   How
%     - Potential of the pattern for 1 A: u = sum_c w_c V_c (out.V_per_A, superposition).
%     - Field in every tetrahedron (elementFieldP1); for each tissue (endoneurium, epineurium,
%       muscle) the elements are sorted by decreasing |E| and their volumes accumulated: the curve
%       "volume of tissue with |E| >= x".
%     - Perineurium: it is a thin layer with no elements (thinLayerInterface); the field inside it is
%       the potential jump across it divided by its real thickness, |u_in - u_out| / t, evaluated at
%       each duplicated node pair.
%     - Mirror: a contact that is not solved but has a solved mirror partner (same type, radius,
%       angle, opposite z; D-12) is replaced by that partner. The |E| distribution of the mirrored
%       pattern is the mirror image, so the exposed volumes are the same.
%     - The contact elements themselves (metal) are left out.
%
%   Inputs
%     out      struct saved by main_pipeline_3d (results/main3d/<id>_3d.mat)
%     ids      contact ids of the pattern
%     weights  fraction of the total current on each contact (sign = polarity); the potential is
%              computed for 1 A of total current
%
%   Output: struct
%     tissue          ["endoneurium" "epineurium" "muscle"]
%     E_perA          1 x 3 cell, |E| [V/m per A] of the elements of each tissue, sorted decreasing
%     cumVol_m3       1 x 3 cell, cumulative volume [m^3] of those elements (same order)
%     periField_perA  field across the perineurium [V/m per A] at each duplicated node, decreasing
%     solvedIds       the solved contacts actually used (after the mirror replacement)
%     mirrored        true where a contact was replaced by its mirror partner

arguments
    out (1,1) struct
    ids (1,:) double
    weights (1,:) double
end
assert(numel(ids) == numel(weights), 'fieldExposure:size', 'ids and weights must have the same length.');

% --- pattern potential for 1 A (mirror replacement where needed) -------------------------------------
sites = out.sites;   solved = double(out.solvedIds(:)');
nAll = size(out.mesh.nodes_m, 2);
u = zeros(nAll, 1);   used = zeros(size(ids));   mirrored = false(size(ids));
for k = 1:numel(ids)
    c = ids(k);
    if ~ismember(c, solved)
        c = mirrorPartner(sites, c);
        assert(ismember(c, solved), 'fieldExposure:notSolved', ...
               'Contact %d is not solved and neither is its mirror partner.', ids(k));
        mirrored(k) = true;
    end
    used(k) = c;
    u = u + weights(k) * double(out.V_per_A(:, solved == c));
end

% --- field in the tissue elements ------------------------------------------------------------------
tets = out.mesh.tets;
elemTissue = out.cells.tissue(double(out.mesh.elemCell));   elemTissue = elemTissue(:)';
[E, vol] = elementFieldP1(out.mesh.nodes_m, tets, u);
Emag = vecnorm(E, 2, 1);
fx = struct();
fx.tissue = ["endoneurium" "epineurium" "muscle"];
codes = [1 3 4];                                              % tissue codes of main_pipeline_3d
fx.E_perA = cell(1, 3);   fx.cumVol_m3 = cell(1, 3);
for k = 1:3
    sel = elemTissue == codes(k);
    [e, order] = sort(Emag(sel), 'descend');
    v = vol(sel);
    fx.E_perA{k} = e;
    fx.cumVol_m3{k} = cumsum(v(order));
end

% --- field across the perineurium (thin layer) --------------------------------------------------------
nOrig = double(out.mesh.nOriginalNodes);
dupOf = double(out.mesh.dupOf(:));
jump = abs(u(nOrig + (1:numel(dupOf))') - u(dupOf));
fx.periField_perA = sort(jump / out.params.perineurium.realThickness_m, 'descend');
fx.solvedIds = used;
fx.mirrored = mirrored;
end


function partner = mirrorPartner(sites, id)
%MIRRORPARTNER Contact with the same type, radius and angle at the opposite z (D-12); NaN if none.
s = sites(sites.id == id, :);
tol = 1e-9;
cand = sites.type == s.type & abs(sites.z_m + s.z_m) < tol & abs(sites.r_m - s.r_m) < tol & ...
       abs(mod(sites.theta_deg - s.theta_deg + 180, 360) - 180) < 1e-6 & sites.id ~= id;
partner = sites.id(cand);
if isempty(partner), partner = NaN; end
end
