function [v, tetId] = interpolateP1(tr, u, points_m)
%INTERPOLATEP1 Value of a P1 (linear-tetrahedra) field at arbitrary points.
%   [v, tetId] = interpolateP1(tr, u, points_m)
%
%   Inputs
%     tr        triangulation(tets', nodes') of the mesh (build it once: it is reused for every
%               field on the same mesh). With the perineurium thin layer, pass the mesh AFTER the
%               node duplication (thinLayerInterface: tl.tets, tl.nodes_m), so that a point inside
%               a fascicle takes the endoneurium-side value.
%     u         nNodes x nFields nodal values
%     points_m  nPoints x 3 query points [m]
%
%   Outputs
%     v      nPoints x nFields interpolated values (NaN outside the mesh)
%     tetId  nPoints x 1 element containing each point (NaN outside)

arguments
    tr triangulation
    u double
    points_m (:,3) double
end
[tetId, bary] = pointLocation(tr, points_m);
v = nan(size(points_m, 1), size(u, 2));
ok = ~isnan(tetId);
nodes = tr.ConnectivityList(tetId(ok), :);            % nOk x 4
for f = 1:size(u, 2)
    uf = u(:, f);
    v(ok, f) = sum(bary(ok, :) .* uf(nodes), 2);
end
end
