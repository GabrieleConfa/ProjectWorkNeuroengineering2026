classdef testLocateRegularGrid < matlab.unittest.TestCase
%TESTLOCATEREGULARGRID locateRegularGrid gives the same interpolation as pointLocation, and the
% block-wise assembly gives the same stiffness matrix as the one-block assembly.
%   Run with:  setup_paths; runtests('tests/testLocateRegularGrid.m')

    methods (Test)
        function sameAsPointLocation(tc)
            % same elements found (up to ties on shared faces) as MATLAB pointLocation: checked through the exact interpolation of a linear field
            rng(3);
            P = [rand(400, 3); dec2bin(0:7) - '0'] * 1e-3;        % random points + cube corners [m]
            T = delaunay(P);   X = P';   T = T';
            xs = linspace(-0.1e-3, 1.1e-3, 13);   ys = linspace(0, 1e-3, 11);   zs = linspace(0.05e-3, 0.95e-3, 7);
            [GX, GY, GZ] = ndgrid(xs, ys, zs);
            mask = GX + GY < 1.5e-3;                              % an arbitrary subset of the grid
            [id1, b1] = locateRegularGrid(X, T, {xs, ys, zs}, mask, 50);
            [id2, b2] = pointLocation(triangulation(T', X'), [GX(mask), GY(mask), GZ(mask)]);
            tc.verifyEqual(isnan(id1), isnan(id2));               % same points outside the mesh
            % a linear field is interpolated exactly by any element containing the point
            f = @(x) 2 + 3*x(1,:) - 5*x(2,:) + 7*x(3,:);
            ok = ~isnan(id1);
            v1 = sum(b1(ok, :) .* reshape(f(X(:, T(:, id1(ok)))), 4, [])', 2);
            v2 = sum(b2(ok, :) .* reshape(f(X(:, T(:, id2(ok)))), 4, [])', 2);
            tc.verifyEqual(v1, v2, 'AbsTol', 1e-12);
            pts = [GX(mask), GY(mask), GZ(mask)]';
            tc.verifyEqual(v1, f(pts(:, ok))', 'AbsTol', 1e-12);
        end

        function blockAssemblySameMatrix(tc)
            % assembling the stiffness matrix in blocks of 97 elements must give the same matrix as one block
            rng(4);
            P = rand(300, 3) * 1e-3;   T = delaunay(P)';   X = P';
            sig = rand(3, size(T, 2));
            [K1, v1] = assembleConductionP1(X, T, sig);
            [K2, v2] = assembleConductionP1(X, T, sig, 97);
            tc.verifyEqual(full(K2), full(K1), 'AbsTol', 1e-12 * max(abs(K1(:))));
            tc.verifyEqual(v2, v1);
        end
    end
end
