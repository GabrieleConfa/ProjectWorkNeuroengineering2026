function [E_Vpm, vol_m3] = elementFieldP1(nodes_m, tets, u_V, blockSize)
%ELEMENTFIELDP1 Electric field E = -grad(V) on each linear tetrahedron of a P1 solution.
%   [E_Vpm, vol_m3] = elementFieldP1(nodes_m, tets, u_V)
%   [E_Vpm, vol_m3] = elementFieldP1(nodes_m, tets, u_V, blockSize)
%
%   Purpose
%     The FEM potential is linear inside each tetrahedron, so its gradient, and the electric field
%     E = -grad V, is constant in each element. This function returns that constant field. It is
%     used by the field-exposure safety check (fieldExposure, scripts/run_field_safety.m).
%
%   How
%     With the barycentric coordinates lambda_i of the element (grad lambda_i constant, computed as
%     in assembleConductionP1): V = sum_i u_i lambda_i, so grad V = sum_i u_i grad lambda_i, and with
%     grad lambda_1 = -(grad lambda_2 + grad lambda_3 + grad lambda_4):
%       grad V = (u_2 - u_1) g_2 + (u_3 - u_1) g_3 + (u_4 - u_1) g_4.
%     Done block by block to keep the temporary arrays small on a 4.5-million-element mesh.
%
%   Caveat: at the sharp edges of a contact the exact field is singular; the P1 field there grows
%   when the mesh is refined. Use volumes above a field level (fieldExposure), not the peak value.
%
%   Inputs
%     nodes_m    3 x nNodes node coordinates [m]
%     tets       4 x nElem node indices of each tetrahedron
%     u_V        nNodes x 1 nodal potential [V]
%     blockSize  optional, elements per block (default 5e5)
%
%   Outputs
%     E_Vpm      3 x nElem electric field [V/m] in each element
%     vol_m3     1 x nElem element volumes [m^3]

arguments
    nodes_m (3,:) double
    tets (4,:) {mustBeInteger, mustBePositive}
    u_V (:,1) double
    blockSize (1,1) double {mustBePositive} = 5e5
end

nE = size(tets, 2);
E_Vpm = zeros(3, nE);   vol_m3 = zeros(1, nE);
for first = 1:blockSize:nE
    blk = first:min(first + blockSize - 1, nE);
    T = double(tets(:, blk));
    p1 = nodes_m(:, T(1, :));
    a = nodes_m(:, T(2, :)) - p1;   b = nodes_m(:, T(3, :)) - p1;   c = nodes_m(:, T(4, :)) - p1;
    V6 = dot(a, cross(b, c));                                 % 6 x signed volume
    g2 = cross(b, c) ./ V6;   g3 = cross(c, a) ./ V6;   g4 = cross(a, b) ./ V6;   % grad lambda_2..4
    u1 = u_V(T(1, :))';   u2 = u_V(T(2, :))';   u3 = u_V(T(3, :))';   u4 = u_V(T(4, :))';
    E_Vpm(:, blk) = -((u2 - u1) .* g2 + (u3 - u1) .* g3 + (u4 - u1) .* g4);
    vol_m3(blk) = abs(V6) / 6;
end
end
