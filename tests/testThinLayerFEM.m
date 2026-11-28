classdef testThinLayerFEM < matlab.unittest.TestCase
%TESTTHINLAYERFEM Analytical checks of the project's P1 assembly and perineurium thin layer.
%   Run with:  setup_paths; runtests('tests/testThinLayerFEM.m')
%   The 3D pipeline (src/main/main_pipeline_3d.m) assembles its own FE system with
%   assembleConductionP1 + thinLayerInterface; these cases have exact solutions that linear elements
%   reproduce to round-off, so any error in the assembly shows up immediately.

    properties (Constant)
        mm = 1e-3
    end

    methods (Test)
        function linearFieldAnisotropic(tc)
            % Block 1 x 1 x 2 mm, sigma = diag(0.1, 0.2, 0.5), V = 0 at z = 0 and 1 V at z = top:
            % exact V = z / L whatever sigma, current = sigma_z * A / L.
            mm = tc.mm;   L = 2*mm;   A = mm^2;   sig = [0.1; 0.2; 0.5];
            g = generateMesh(fegeometry(multicuboid(mm, mm, L)), Hmax=0.4*mm, GeometricOrder="linear");
            X = g.Mesh.Nodes;   T = g.Mesh.Elements;
            K = assembleConductionP1(X, T, sig);
            bot = find(abs(X(3,:)) < 1e-12);   top = find(abs(X(3,:) - L) < 1e-12);
            u = solveDirichlet(K, zeros(size(X,2),1), [bot top], [zeros(size(bot)) ones(size(top))]);
            tc.verifyLessThan(max(abs(u' - X(3,:) / L)), 1e-10);
            Itop = sum(K(top, :) * u);                % current injected at the top electrode
            tc.verifyEqual(Itop, sig(3) * A / L, 'RelTol', 1e-9);
        end

        function thinLayerJump(tc)
            % Two stacked blocks, inside (bottom, h1, sigma1) and outside (top, h2, sigma2), with a
            % thin layer of sheet conductance Gs at z = h1. V = 0 at the bottom, 1 V at the top.
            % Series resistances per unit area: J = 1 / (h1/s1 + 1/Gs + h2/s2), jump = J / Gs.
            mm = tc.mm;   h1 = 1*mm;   h2 = 1.5*mm;   s1 = 0.3;   s2 = 0.1;   Gs = 63;
            g = fegeometry(multicuboid(mm, mm, [h1 h2], ZOffset=[0 h1]));
            g = generateMesh(g, Hmax=0.35*mm, GeometricOrder="linear");
            X = g.Mesh.Nodes;   T = g.Mesh.Elements;   nE = size(T, 2);
            zc = mean(reshape(X(3, T), 4, nE), 1);
            inside = zc < h1;
            iface = find(abs(X(3,:) - h1) < 1e-12);
            tl = thinLayerInterface(X, T, inside, iface, Gs);
            sig = repmat([s2; s2; s2], 1, nE);   sig(:, inside) = s1;
            K = assembleConductionP1(tl.nodes_m, tl.tets, sig) + tl.K;
            Xa = tl.nodes_m;
            bot = find(abs(Xa(3,:)) < 1e-12);   top = find(abs(Xa(3,:) - (h1 + h2)) < 1e-12);
            u = solveDirichlet(K, zeros(tl.nNodes,1), [bot top], [zeros(size(bot)) ones(size(top))]);
            J = 1 / (h1/s1 + 1/Gs + h2/s2);
            tc.verifyEqual(u(tl.copyOf(iface)), repmat(J*h1/s1, numel(iface), 1), 'AbsTol', 1e-9);   % inside
            tc.verifyEqual(u(iface), repmat(1 - J*h2/s2, numel(iface), 1), 'AbsTol', 1e-9);        % outside
            tc.verifyEqual(sum(tl.area_m2), mm^2, 'RelTol', 1e-9);       % the whole interface is found
            tc.verifyEqual(sum(K(top,:) * u), J * mm^2, 'RelTol', 1e-8);  % current
        end

        function volumeSourceConservation(tc)
            % A current I injected as a uniform source in an inner cube leaves entirely through the
            % grounded outer faces.
            mm = tc.mm;   I = 1e-6;
            g = fegeometry(multicuboid(4*mm, 4*mm, 4*mm));
            g = addCell(g, translate(fegeometry(multicuboid(mm, mm, mm)), [0 0 1.5*mm]));
            g = generateMesh(g, Hmax=0.5*mm, GeometricOrder="linear");
            X = g.Mesh.Nodes;   T = g.Mesh.Elements;
            inner = findCell(g, [0 0 2*mm]);
            src = findElements(g.Mesh, "region", "Cell", inner);
            c = mean(reshape(X(1, T(:, src)), 4, []), 1);
            tc.assertLessThan(max(abs(c)), 0.6*mm);  % the source cell is the inner cube
            [K, vol] = assembleConductionP1(X, T, [0.2; 0.2; 0.5]);
            F = accumarray(reshape(T(:, src), [], 1), reshape(repmat(I / sum(vol(src)) * vol(src) / 4, 4, 1), [], 1), ...
                           [size(X,2) 1]);
            outer = find(abs(abs(X(1,:)) - 2*mm) < 1e-12 | abs(abs(X(2,:)) - 2*mm) < 1e-12 | ...
                         abs(X(3,:)) < 1e-12 | abs(X(3,:) - 4*mm) < 1e-12);
            u = solveDirichlet(K, F, outer, zeros(size(outer)));
            Iout = -sum(K(outer, :) * u - F(outer));
            tc.verifyEqual(Iout, I, 'RelTol', 1e-9);
            tc.verifyEqual(sum(F), I, 'RelTol', 1e-12);
            tc.verifyGreaterThan(min(u), -1e-15);     % source only: V >= 0 everywhere
        end

        function interpolationExactForLinear(tc)
            % interpolateP1 must reproduce a linear field exactly at arbitrary points
            mm = tc.mm;
            g = generateMesh(fegeometry(multicuboid(mm, mm, mm)), Hmax=0.4*mm, GeometricOrder="linear");
            X = g.Mesh.Nodes;   T = g.Mesh.Elements;
            u = (2*X(1,:) - 3*X(2,:) + 5*X(3,:))';
            tr = triangulation(double(T'), X');
            P = [0.1 0.2 0.3; -0.3 0.4 0.9; 0.45 -0.45 0.05] * mm;
            v = interpolateP1(tr, u, P);
            tc.verifyEqual(v, 2*P(:,1) - 3*P(:,2) + 5*P(:,3), 'AbsTol', 1e-15);
            tc.verifyTrue(isnan(interpolateP1(tr, u, [5 5 5]*mm)));
        end
    end
end

function u = solveDirichlet(K, F, dNodes, dVal)
%SOLVEDIRICHLET Solve K u = F with u(dNodes) = dVal.
n = size(K, 1);   u = zeros(n, 1);   u(dNodes) = dVal;
free = setdiff(1:n, dNodes);
u(free) = K(free, free) \ (F(free) - K(free, dNodes) * u(dNodes));
end
