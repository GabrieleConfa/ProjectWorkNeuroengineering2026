classdef testRandomFascicleShapes < matlab.unittest.TestCase % Checks of the random irregular fascicle generator.
    methods (Test)
        function areaInsideAndGapKept(tc) % Exact area, inside the nerve, no overlap, reproducible.
            q = quicklook3DParams(); ip = insertionParams(); % Same inputs as the study.
            S = randomFascicleShapes(q,ip,ip.seed); n = numel(S.polygons); % Generate the shapes.
            for k = 1:n, tc.verifyEqual(polyarea(S.polygons{k}(:,1),S.polygons{k}(:,2)),S.area_m2(k),'RelTol',1e-6); end % Area is preserved.
            inner = polybuffer(polyshape(S.nerve),-q.geometry.tPeri_m); % Region inside the perineurium margin.
            for k = 1:n, tc.verifyTrue(all(inpolygon(S.polygons{k}(:,1),S.polygons{k}(:,2),inner.Vertices(:,1),inner.Vertices(:,2)))); end % Inside the nerve.
            for i = 1:n-1, for j = i+1:n, tc.verifyFalse(overlaps(polyshape(S.polygons{i}),polyshape(S.polygons{j}))); end, end % No overlap.
            T = randomFascicleShapes(q,ip,ip.seed); tc.verifyEqual(T.roundness,S.roundness); % Same seed, same shapes.
            tc.verifyGreaterThan(mean(S.roundness),0.55); tc.verifyLessThan(mean(S.roundness),0.85); % Statistics near the histology value.
        end
    end
end
