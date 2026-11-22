function [tetId, bary] = locateRegularGrid(X, T, axes_m, mask, batchSize)
%LOCATEREGULARGRID Find the tetrahedron containing each point of a regular 3D grid.
%   [tetId, bary] = locateRegularGrid(X, T, axes_m, mask)
%   [tetId, bary] = locateRegularGrid(X, T, axes_m, mask, batchSize)
%
%   Same result as pointLocation(triangulation(T', X'), P) for the grid points P selected by mask,
%   but much faster for a regular grid: instead of searching the mesh for every point, every
%   tetrahedron lists the few grid points inside its bounding box (index arithmetic, since the grid
%   is regular) and keeps those whose barycentric coordinates are all >= 0. A point on a face shared
%   by two elements gets the first one found (any of them interpolates a P1 field correctly).
%   Adapted from Marco's branch (cuff-oval-viewer-convergence, 2026-11-22); checked against
%   pointLocation in tests/testLocateRegularGrid.m.
%
%   Inputs
%     X          3 x nNodes node coordinates [m]
%     T          4 x nElem tetrahedra (pass only the elements that can contain the points, e.g. the
%                nerve and contact elements, to save time)
%     axes_m     {xs, ys, zs}: the uniformly spaced grid coordinates along x, y, z [m]
%     mask       logical array numel(xs) x numel(ys) x numel(zs): the points to locate
%     batchSize  optional, elements processed at a time (memory), default 20000
%   Outputs
%     tetId      (nnz(mask),1) column of T containing each point, NaN if none; order = find(mask)
%     bary       (nnz(mask),4) barycentric coordinates in that element (NaN if none)

arguments
    X (3,:) double
    T (4,:) double
    axes_m (1,3) cell
    mask logical
    batchSize (1,1) double {mustBePositive, mustBeInteger} = 20000
end

dims = cellfun(@numel, axes_m);
origin = cellfun(@(v) v(1), axes_m);
step = cellfun(@(v) v(2) - v(1), axes_m);
assert(isequal(size(mask, 1:3), dims), 'locateRegularGrid:mask', 'mask must be numel(xs) x numel(ys) x numel(zs).');
ids = nan(prod(dims), 1);   weights = nan(prod(dims), 4);
tol = 1e-10;
for first = 1:batchSize:size(T, 2)
    el = first:min(first + batchSize - 1, size(T, 2));
    corners = T(:, el);
    lo = zeros(numel(el), 3);   hi = lo;                     % grid index range of each bounding box
    for d = 1:3
        coord = reshape(X(d, corners), 4, []);
        lo(:, d) = max(1, ceil((min(coord, [], 1)' - origin(d)) / step(d) - tol) + 1);
        hi(:, d) = min(dims(d), floor((max(coord, [], 1)' - origin(d)) / step(d) + tol) + 1);
    end
    span = max(0, hi - lo + 1);   counts = prod(span, 2);
    if ~any(counts), continue; end
    row = repelem((1:numel(el))', counts);                  % candidate -> element of the batch
    offset = (0:numel(row)-1)' - repelem(cumsum(counts) - counts, counts);
    ijk = lo(row, :) + [mod(offset, span(row,1)), mod(floor(offset ./ span(row,1)), span(row,2)), ...
                        floor(offset ./ (span(row,1) .* span(row,2)))];
    k = sub2ind(dims, ijk(:,1), ijk(:,2), ijk(:,3));
    keep = mask(k) & isnan(ids(k));                          % wanted and not found yet
    k = k(keep);   row = row(keep);   ijk = ijk(keep, :);
    if isempty(k), continue; end
    pts = origin + (ijk - 1) .* step;
    tc = corners(:, row);   p1 = X(:, tc(1,:));
    a = X(:, tc(2,:)) - p1;   b = X(:, tc(3,:)) - p1;   c = X(:, tc(4,:)) - p1;
    dp = pts' - p1;   V6 = dot(a, cross(b, c, 1), 1);
    w = [dot(dp, cross(b, c, 1), 1); dot(dp, cross(c, a, 1), 1); dot(dp, cross(a, b, 1), 1)] ./ V6;
    w = [1 - sum(w, 1); w]';
    inside = all(w >= -tol & w <= 1 + tol, 2);
    % a point can be in two elements of the same batch: keep the first
    [kIn, firstOcc] = unique(k(inside), 'stable');
    rIn = row(inside);   wIn = w(inside, :);
    ids(kIn) = el(rIn(firstOcc));   weights(kIn, :) = wIn(firstOcc, :);
end
tetId = ids(mask);   bary = weights(mask, :);
end
