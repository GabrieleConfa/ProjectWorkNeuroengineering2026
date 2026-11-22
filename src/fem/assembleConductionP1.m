function [K, vol_m3] = assembleConductionP1(nodes_m, tets, sigma_Spm)
%ASSEMBLECONDUCTIONP1 Stiffness matrix of div(sigma grad V) = 0 on linear tetrahedra.
%   [K, vol_m3] = assembleConductionP1(nodes_m, tets, sigma_Spm)
%
%   Standard P1 finite elements: on each tetrahedron V is linear, grad(lambda_i) is constant and
%   K_e(i,j) = vol * grad(lambda_i)' * S * grad(lambda_j), with S = diag(sigma_x, sigma_y, sigma_z)
%   (diagonal conductivity tensor: z = nerve axis = fibre direction).
%   Written for the project because the perineurium thin layer (thinLayerInterface) needs the nodes
%   of the fascicle boundary to be split, which the PDE Toolbox assembly does not allow.
%
%   Inputs
%     nodes_m    3 x nNodes node coordinates [m]
%     tets       4 x nElem node indices of each linear tetrahedron
%     sigma_Spm  3 x nElem (or 3 x 1, same for all) conductivities [S/m]: sigma_x, sigma_y, sigma_z
%
%   Outputs
%     K       nNodes x nNodes sparse stiffness matrix [S]: K*V = injected currents [A]
%     vol_m3  1 x nElem element volumes [m^3]

arguments
    nodes_m (3,:) double
    tets (4,:) {mustBeInteger, mustBePositive}
    sigma_Spm (3,:) double
end

nN = size(nodes_m, 2);   nE = size(tets, 2);
if size(sigma_Spm, 2) == 1, sigma_Spm = repmat(sigma_Spm, 1, nE); end
assert(size(sigma_Spm, 2) == nE, 'assembleConductionP1:sigma', 'sigma must be 3 x 1 or 3 x nElem.');

% Bound transient memory for a whole-nerve mesh on an 8 GB machine.
blockSize = 150000; % Limit the 16 local stiffness entries per tetrahedron.
if nE > blockSize
    K = sparse(nN,nN); vol_m3 = zeros(1,nE); % Accumulate into the common global node numbering.
    for first = 1:blockSize:nE
        block = first:min(first+blockSize-1,nE); % Select the next contiguous element block.
        [part,vol_m3(block)] = assembleConductionP1(nodes_m,tets(:,block),sigma_Spm(:,block)); % Assemble the same element formula.
        K = K + part; % Sum contributions at shared nodes across blocks.
    end
    return
end
p1 = nodes_m(:, tets(1,:));
a = nodes_m(:, tets(2,:)) - p1;   b = nodes_m(:, tets(3,:)) - p1;   c = nodes_m(:, tets(4,:)) - p1;
V6 = dot(a, cross(b, c));                             % 6 x signed volume
assert(all(abs(V6) > 0), 'assembleConductionP1:degenerate', 'Degenerate (zero-volume) element.');
g = zeros(3, 4, nE);                                  % grad(lambda_i), i = 1..4
g(:, 2, :) = cross(b, c) ./ V6;
g(:, 3, :) = cross(c, a) ./ V6;
g(:, 4, :) = cross(a, b) ./ V6;
g(:, 1, :) = -(g(:, 2, :) + g(:, 3, :) + g(:, 4, :));
vol_m3 = abs(V6) / 6;

I = zeros(16, nE);   J = zeros(16, nE);   Vals = zeros(16, nE);
k = 0;
for i = 1:4
    gi = reshape(g(:, i, :), 3, nE) .* sigma_Spm;     % S * grad(lambda_i)
    for j = 1:4
        k = k + 1;
        gj = reshape(g(:, j, :), 3, nE);
        I(k, :) = tets(i, :);   J(k, :) = tets(j, :);
        Vals(k, :) = vol_m3 .* sum(gi .* gj, 1);
    end
end
K = sparse(I(:), J(:), Vals(:), nN, nN);
end
