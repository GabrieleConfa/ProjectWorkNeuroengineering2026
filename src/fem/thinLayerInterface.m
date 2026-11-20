function tl = thinLayerInterface(nodes_m, tets, insideElem, interfaceNodes, Gs_Spm2)
%THINLAYERINTERFACE Perineurium as a thin resistive layer between two regions of a P1 mesh.
%   tl = thinLayerInterface(nodes_m, tets, insideElem, interfaceNodes, Gs_Spm2)
%
%   Physics: the perineurium is too thin to mesh (14 um along 70 mm), but electrically it only acts
%   as a sheet resistance rho*t: the current density crossing it is J = Gs (V_in - V_out), with
%   Gs = 1 / (rho*t) = sigma_peri / t_peri [S/m^2], and the potential JUMPS across it. This is the
%   "contact impedance" / thin-layer approximation (Pelot 2019; ASCENT, Davis 2023).
%
%   Implementation: every node on the interface is duplicated. The elements on the inside
%   (endoneurium) use the copy, the other elements keep the original node, so the two sides are
%   disconnected; they are then coupled through the interface "membrane" matrix
%     K_int = Gs * [M -M; -M M],   M = surface mass matrix of the interface triangles
%   (for a P1 triangle of area A: M = A/12 * [2 1 1; 1 2 1; 1 1 2]).
%
%   Inputs
%     nodes_m         3 x nNodes node coordinates [m]
%     tets            4 x nElem linear tetrahedra
%     insideElem      1 x nElem logical, true for the elements inside the layer (endoneurium)
%     interfaceNodes  node indices on the layer (nodes shared with a contact should be left out:
%                     where a contact pierces the layer the two sides are shorted by the metal)
%     Gs_Spm2         sheet conductance of the layer [S/m^2]
%
%   Output: struct
%     tets      4 x nElem tetrahedra with the inside elements renumbered to the copies
%     nNodes    number of nodes after the duplication (original + copies)
%     nodes_m   3 x nNodes coordinates (the copies have the coordinates of their originals)
%     dupOf     (nCopies x 1) original node of each copy; copy k has index nOrig + k
%     copyOf    (nOrig x 1) index of the copy of each original node, 0 if not duplicated
%     K         nNodes x nNodes sparse interface matrix [S]
%     tri       3 x nTri interface triangles (original, outside numbering)
%     area_m2   1 x nTri their areas

arguments
    nodes_m (3,:) double
    tets (4,:) double
    insideElem (1,:) logical
    interfaceNodes (:,1) double
    Gs_Spm2 (1,1) double {mustBePositive}
end

nN = size(nodes_m, 2);
interfaceNodes = unique(interfaceNodes);
nC = numel(interfaceNodes);
copyOf = zeros(nN, 1);   copyOf(interfaceNodes) = nN + (1:nC)';
onLayer = false(nN, 1);  onLayer(interfaceNodes) = true;

% Interface triangles: faces of inside elements that belong to only one inside element (boundary of
% the inside region) and whose three nodes are all on the layer.
T = tets(:, insideElem);
F = [T([1 2 3], :), T([1 2 4], :), T([1 3 4], :), T([2 3 4], :)];
Fs = sort(F, 1)';
[uF, ~, ic] = unique(Fs, 'rows');
cnt = accumarray(ic, 1);
tri = uF(cnt == 1 & all(onLayer(uF), 2), :)';

e1 = nodes_m(:, tri(2,:)) - nodes_m(:, tri(1,:));
e2 = nodes_m(:, tri(3,:)) - nodes_m(:, tri(1,:));
area_m2 = 0.5 * vecnorm(cross(e1, e2), 2, 1);

nAll = nN + nC;
Mloc = [2 1 1; 1 2 1; 1 1 2] / 12;
I = zeros(36, size(tri, 2));   J = I;   V = I;
k = 0;
for a = 1:3
    for b = 1:3
        m = Gs_Spm2 * area_m2 * Mloc(a, b);
        o_a = tri(a, :);   o_b = tri(b, :);
        i_a = copyOf(o_a)'; i_b = copyOf(o_b)';
        k = k + 1;  I(k,:) = o_a;  J(k,:) = o_b;  V(k,:) =  m;   % out-out
        k = k + 1;  I(k,:) = i_a;  J(k,:) = i_b;  V(k,:) =  m;   % in-in
        k = k + 1;  I(k,:) = o_a;  J(k,:) = i_b;  V(k,:) = -m;   % out-in
        k = k + 1;  I(k,:) = i_a;  J(k,:) = o_b;  V(k,:) = -m;   % in-out
    end
end
tl.K = sparse(I(:), J(:), V(:), nAll, nAll);

tetsNew = tets;
Ti = tetsNew(:, insideElem);
dupl = onLayer(Ti);
Ti(dupl) = copyOf(Ti(dupl));
tetsNew(:, insideElem) = Ti;

tl.tets = tetsNew;
tl.nNodes = nAll;
tl.nodes_m = [nodes_m, nodes_m(:, interfaceNodes)];
tl.dupOf = interfaceNodes;
tl.copyOf = copyOf;
tl.tri = tri;
tl.area_m2 = area_m2;
end
