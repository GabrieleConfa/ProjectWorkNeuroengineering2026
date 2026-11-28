classdef testFieldSafety < matlab.unittest.TestCase
%TESTFIELDSAFETY Known answers for the field-exposure safety check (src/safety/).
%   - elementFieldP1: a linear potential V = c0 + g.x has the exact field E = -g in every element,
%     and the element volumes add up to the volume of the meshed region;
%   - exposedVolume: reads the cumulative curve correctly and scales with the current (linearity).
%   Run with:  setup_paths; runtests('tests/testFieldSafety.m')

    methods (Test)
        function linearPotentialExactField(tc)
            % a linear potential has the exact field -grad V in every element; the volumes add up to the cube
            rng(5);
            P = [rand(300, 3); dec2bin(0:7) - '0'] * 1e-3;     % points in a 1 mm cube + its corners [m]
            T = delaunay(P)';   X = P';
            g = [120; -45; 300];                              % potential gradient [V/m]
            u = 0.7 + g' * X;                                 % V = c0 + g.x at the nodes [V]
            [E, vol] = elementFieldP1(X, T, u(:), 200);       % small blocks: exercises the block loop
            tc.verifyEqual(E, repmat(-g, 1, size(T, 2)), 'AbsTol', 1e-8 * norm(g));
            tc.verifyEqual(sum(vol), 1e-9, 'RelTol', 1e-10); % the cube volume [m^3]
        end

        function exposedVolumeReadsCurve(tc)
            % exposed volume read from the cumulative curve, and scaled with the current (linear model)
            fx = struct('tissue', "t", 'E_perA', {{[30 20 10]}}, 'cumVol_m3', {{[1 3 6]}});
            tc.verifyEqual(exposedVolume(fx, 20, 1), 3);     % two elements have |E| >= 20 V/m per A
            tc.verifyEqual(exposedVolume(fx, 20, 2), 6);     % at 2 A all three exceed 20 V/m
            tc.verifyEqual(exposedVolume(fx, 100, 1), 0);    % none
            tc.verifyEqual(exposedVolume(fx, 20, [1; 2]), [3; 6]);
        end
    end
end
