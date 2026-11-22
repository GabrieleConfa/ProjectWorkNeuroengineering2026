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
            g = tc.p.geometry.generator;   f = tc.P1.fascicles;   ax = tc.P1.nerveSemiAxes_m;
            tc.verifyEqual(height(f), g.nFascicles);
            tc.verifyEqual(pi * tc.P1.nerveRadius_m^2, g.nerveArea_m2, 'RelTol', 1e-12);
            tc.verifyEqual(pi * ax(1) * ax(2), g.nerveArea_m2, 'RelTol', 1e-12);     % ellipse area
            if g.nerveShape == "ellipse"
                tc.verifyEqual(ax(2) / ax(1), g.nerveRoundness, 'RelTol', 1e-12);
            end
            tc.verifyEqual(sum(pi * f.rEndo_m.^2), g.totalFascicleArea_m2, 'RelTol', 1e-12);
            tc.verifyEqual(pi * f.rEndo_m(1)^2, g.largestFascicleArea_m2, 'RelTol', 1e-12);
            tc.verifyTrue(all(f.rEndo_m(2:end) < f.rEndo_m(1)));
        end

        function noOverlapAndInside(tc)
            % inside: distance of every fascicle centre from the nerve outline (sampled much more
            % finely than in the generator) >= outer radius + gap, up to the sampling error (< 0.05 um)
            g = tc.p.geometry.generator;   f = tc.P1.fascicles;
            rOut = f.rEndo_m + f.tPeri_m;
            bnd = nerveOutline(tc.P1.nerveSemiAxes_m, tc.P1.nerveRotation_deg, 40000);
            for i = 1:height(f)
                c = [f.x_m(i), f.y_m(i)];
                tc.verifyLessThan(ellipseQ(c, tc.P1.nerveSemiAxes_m, tc.P1.nerveRotation_deg), 1);
                tc.verifyGreaterThanOrEqual(min(vecnorm(bnd - c, 2, 2)), rOut(i) + g.minGapNerve_m - 0.05e-6);
            end
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

        function circleOptionStillWorks(tc)
            % nerveShape = "circle" gives back the round nerve of the first 3D model (D-09)
            q = tc.p;   q.geometry.generator.nerveShape = "circle";
            Pc = makePatient("TC", 1, q);
            tc.verifyEqual(Pc.nerveSemiAxes_m, Pc.nerveRadius_m * [1 1]);
            rOut = Pc.fascicles.rEndo_m + Pc.fascicles.tPeri_m;
            tc.verifyLessThanOrEqual(max(hypot(Pc.fascicles.x_m, Pc.fascicles.y_m) + rOut), ...
                                     Pc.nerveRadius_m - q.geometry.generator.minGapNerve_m + 1e-12);
        end
    end
end
