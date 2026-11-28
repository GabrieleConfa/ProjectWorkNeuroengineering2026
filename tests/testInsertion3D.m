classdef testInsertion3D < matlab.unittest.TestCase
%TESTINSERTION3D The lead-insertion model (Marco's simulateInsertion) on the patient of the main 3D
% model, with the two pierce models.
%   - reproducible: the same random stream gives the same lead path;
%   - a completed insertion ends on the far side of the nerve, against the wall;
%   - probability model, limits of the law: P -> 1 pierces every fascicle touched, P -> 0 none;
%   - the law itself: grows with the radius, falls with the angle, the values quoted in its docstring.
%   Run with:  setup_paths; runtests('tests/testInsertion3D.m')

    properties
        q; ip; polys; st;
    end

    methods (TestClassSetup)
        function setup(tc)
            % patient of the main 3D model, its circular fascicle outlines, the nerve outline and the 8 stations
            p = pw6_3d_params();
            patient = makePatient(p.patient.id, p.patient.seed, p);
            tc.q = insertionGeometry(patient, p);
            tc.ip = insertionParams();
            th = (0:127)' * 2 * pi / 128;   F = tc.q.geometry.fascicles_m;
            tc.polys = arrayfun(@(k) F(k, 1:2) + F(k, 3) * [cos(th) sin(th)], (1:size(F, 1))', 'UniformOutput', false);
            nerve = nerveOutline(patient.nerveSemiAxes_m, patient.nerveRotation_deg, 720);
            tc.st = insertionStations(tc.q, nerve);
        end
    end

    methods (Test)
        function reproducible(tc)
            % the same random stream must give exactly the same lead path
            o1 = simulateInsertion(tc.q, tc.ip, tc.polys, tc.st{1}, RandStream('mt19937ar', 'Seed', 7));
            o2 = simulateInsertion(tc.q, tc.ip, tc.polys, tc.st{1}, RandStream('mt19937ar', 'Seed', 7));
            tc.verifyEqual(o1.path, o2.path);
        end

        function completedLeadEndsOnFarWall(tc)
            % Most random insertions do NOT complete (the tip stops against a fascicle next to the entry
            % wall, or gets trapped between two fascicles it does not pierce): insertionScenario rejects
            % them and draws again, with one stream per attempt. Same here: draw until 3 insertions have
            % completed (with the default parameters about 1 attempt in 20-40 completes on UA01), then
            % check that each one crosses the nerve from the entry wall to the far wall.
            % (Until 2026-11-28 this test drew 20 trials from one stream and none completed.)
            nDone = 0;   u = tc.st{1}.u;
            for attempt = 1:400
                rs = RandStream('mt19937ar', 'Seed', tc.ip.seed + 1000 + attempt);
                o = simulateInsertion(tc.q, tc.ip, tc.polys, tc.st{1}, rs);
                if ~o.completed, continue; end
                nDone = nDone + 1;
                tc.verifyLessThan(o.path(1, :) * u', 0);       % enters on the -u side
                tc.verifyGreaterThan(o.path(end, :) * u', 0);  % ends on the +u side
                if nDone == 3, break; end
            end
            tc.verifyEqual(nDone, 3, 'fewer than 3 completed insertions in 400 attempts');
        end

        function probabilityModelLimits(tc)
            % probability pierce model at its limits: r50 -> 0 pierces every fascicle touched, r50 -> large pierces none (and never displaces)
            ipR = tc.ip;   ipR.pierceModel = "probability";   ipR.entryAngleRange_deg = 0;   ipR.entryOffsetSD_m = 0;
            ipR.pierceAngleExponent = 0;                       % no angle dependence here: P = P_r(r)
            ipAll = ipR;   ipAll.pierceR50_m = 1e-9;            % P(pierce) -> 1
            ipNone = ipR;  ipNone.pierceR50_m = 1;              % P(pierce) -> 0
            nAll = 0;
            for iS = 1:numel(tc.st)
                oA = simulateInsertion(tc.q, ipAll, tc.polys, tc.st{iS}, RandStream('mt19937ar', 'Seed', iS));
                tc.verifyEqual(oA.pierce, oA.pierce | oA.deflect); % whatever is touched is pierced
                nAll = nAll + nnz(oA.pierce);
                oN = simulateInsertion(tc.q, ipNone, tc.polys, tc.st{iS}, RandStream('mt19937ar', 'Seed', iS));
                tc.verifyFalse(any(oN.pierce));
                tc.verifyFalse(any(oN.displace));              % this model never displaces
            end
            tc.verifyGreaterThan(nAll, 0);                      % the section is full of fascicles
        end

        function lawShape(tc)
            % pierceProbability: the values quoted in its header, growing with the radius, falling with the angle, halved at 45 deg
            ipL = insertionParams();
            r = [0.173 0.278 0.434 0.689] * 1e-3;              % F10, F7, F3, F1 of UA01
            tc.verifyEqual(pierceProbability(r, 0, ipL), [0.0088 0.133 0.689 0.973], 'AbsTol', 2e-3);
            tc.verifyTrue(all(diff(pierceProbability(r, 0, ipL)) > 0));          % grows with the radius
            g = [0 30 45 60 89];
            tc.verifyTrue(all(diff(pierceProbability(0.5e-3, g, ipL)) < 0));     % falls with the angle
            tc.verifyEqual(pierceProbability(0.5e-3, 45, ipL), 0.5 * pierceProbability(0.5e-3, 0, ipL), 'RelTol', 1e-12);
        end
    end
end
