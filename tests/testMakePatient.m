classdef testMakePatient < matlab.unittest.TestCase
%TESTMAKEPATIENT Checks on the upper-arm cross-section generator (makePatient / generateFascicleLayout).
%   Run with:  setup_paths; runtests('tests/testMakePatient.m')

    properties
        p
        P1
    end

    methods (TestClassSetup)
        function make(tc)
            tc.p = pw6_3d_params();
            tc.P1 = makePatient("T1", 1, tc.p);
        end
    end

    methods (Test)
        function reproducibleFromSeed(tc)
            P1b = makePatient("T1", 1, tc.p);
            tc.verifyEqual(P1b.fascicles, tc.P1.fascicles);
        end

        function literatureTotals(tc)
            g = tc.p.geometry.generator;   f = tc.P1.fascicles;
            tc.verifyEqual(height(f), g.nFascicles);
            tc.verifyEqual(pi * tc.P1.nerveRadius_m^2, g.nerveArea_m2, 'RelTol', 1e-12);
            tc.verifyEqual(sum(pi * f.rEndo_m.^2), g.totalFascicleArea_m2, 'RelTol', 1e-12);
            tc.verifyEqual(pi * f.rEndo_m(1)^2, g.largestFascicleArea_m2, 'RelTol', 1e-12);
            tc.verifyTrue(all(f.rEndo_m(2:end) < f.rEndo_m(1)));
        end

        function noOverlapAndInside(tc)
            g = tc.p.geometry.generator;   f = tc.P1.fascicles;   R = tc.P1.nerveRadius_m;
            rOut = f.rEndo_m + f.tPeri_m;
            tc.verifyLessThanOrEqual(max(hypot(f.x_m, f.y_m) + rOut), R - g.minGapNerve_m + 1e-12);
            for i = 1:height(f)
                for j = i+1:height(f)
                    tc.verifyGreaterThanOrEqual(hypot(f.x_m(i) - f.x_m(j), f.y_m(i) - f.y_m(j)), ...
                        rOut(i) + rOut(j) + g.minGapFascicles_m - 1e-12);
                end
            end
        end

        function notSymmetric(tc)
            % the layout must not be invariant under the mirror x -> -x or y -> -y
            f = tc.P1.fascicles;
            for M = {[-1 1], [1 -1], [-1 -1]}
                m = M{1};
                d = min(hypot(f.x_m - m(1)*f.x_m', f.y_m - m(2)*f.y_m'), [], 2);
                tc.verifyGreaterThan(max(d), 0.1e-3);
            end
        end

        function differentSeedsDiffer(tc)
            P2 = makePatient("T2", 2, tc.p);
            tc.verifyNotEqual(P2.fascicles.x_m, tc.P1.fascicles.x_m);
        end
    end
end
