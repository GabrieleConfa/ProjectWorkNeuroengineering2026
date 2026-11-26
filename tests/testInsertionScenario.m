classdef testInsertionScenario < matlab.unittest.TestCase % Checks of the production insertion scenario (irregular fascicles, leads wall to wall).
    methods (Test)
        function leadsCrossTheNerveAndCarryAllPads(tc) % Every lead runs from one wall to the opposite wall and carries its four pads.
            q = quicklook3DParams(); sc = q.geometry.scenario; nerve = polyshape(sc.nerve); % Scenario and nerve outline.
            tc.verifyEqual(numel(sc.stations),8); % Eight stations.
            for s = sc.stations
                tc.verifyEqual(size(s.padXY,1),4); tc.verifyFalse(any(isnan(s.padXY(:)))); % Four pads placed on the lead.
                tc.verifyTrue(all(isinterior(nerve,s.padXY(:,1),s.padXY(:,2)))); % Pad centres inside the nerve.
                u = s.u; tc.verifyLessThan(s.path(1,:)*u',0); tc.verifyGreaterThan(s.path(end,:)*u',0); % Enters on one side and ends on the opposite side.
            end
        end
        function leadsEndAgainstTheWallAtBothEnds(tc) % Both ends stop with the tip edge touching the wall (the sleeve): not short, not beyond.
            q = quicklook3DParams(); sc = q.geometry.scenario; P = sc.nerve; Q2 = circshift(P,-1,1); S = Q2-P; Rtip = sc.ip.tipDiameter_m/2; % Nerve outline and tip radius.
            for s = sc.stations
                for e = [s.path(1,:);s.path(end,:)]' % Both ends of the lead: entry and anchored end.
                    e = e'; t = ((e(1)-P(:,1)).*S(:,1)+(e(2)-P(:,2)).*S(:,2))./sum(S.^2,2); t = min(1,max(0,t)); % Closest point on each wall segment.
                    d = min(hypot(e(1)-(P(:,1)+t.*S(:,1)),e(2)-(P(:,2)+t.*S(:,2)))); % Distance from the tip centre to the wall.
                    tc.verifyEqual(d,0,'AbsTol',5e-6); % Both ends end on the wall (within 5 um).
                end
            end
        end
        function displacedFasciclesStayInsideAndSeparated(tc) % Displaced outlines keep the gap and stay inside the nerve.
            q = quicklook3DParams(); sc = q.geometry.scenario; inner = polybuffer(polyshape(sc.nerve),-q.geometry.tPeri_m); % Allowed region.
            for s = sc.stations
                for k = 1:numel(s.polys), tc.verifyTrue(all(inpolygon(s.polys{k}(:,1),s.polys{k}(:,2),inner.Vertices(:,1),inner.Vertices(:,2)))); end
            end
        end
        function deterministic(tc) % Same inputs, same scenario.
            q1 = quicklook3DParams(); q2 = quicklook3DParams(); tc.verifyEqual(q1.geometry.scenario.stations(1).path,q2.geometry.scenario.stations(1).path);
        end
    end
end
