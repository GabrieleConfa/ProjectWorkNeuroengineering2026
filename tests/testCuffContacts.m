classdef testCuffContacts < matlab.unittest.TestCase
%TESTCUFFCONTACTS Rigid flat cuff pads on the elliptical nerve (placeCuffContacts, D-11).
%   Run with:  setup_paths; runtests('tests/testCuffContacts.m')
%   Idea of the checks (after Marco's testCuffTangent, 2026-11-22): the face of every pad must be
%   tangent to the native outline at the pad's landmark, both side edges must be exactly edgeIndent
%   below the native surface (the whole face touches the nerve), and on a circle the old radial
%   placement must come back.

    properties
        sites
        edge = 7e-6
    end

    methods (TestClassSetup)
        function load(tc)
            % the provided 48-contact layout, with the 50 um contact thickness used by the 3D model
            p = pw6_3d_params();
            tc.sites = loadActiveSites(p.electrode.sitesFile, 1.7e-3);
            tc.sites.thickness_m(:) = p.electrode.thicknessOverride_m;
        end
    end

    methods (Test)
        function circleGivesRadialPlacement(tc)
            % on a circle the tangent pad must sit radially, as in the original layout: centre depth = edge indent + sagitta of the pad chord
            R = 1.711e-3;
            pat = struct('nerveSemiAxes_m', [R R], 'nerveRotation_deg', 0);
            [s, fit] = placeCuffContacts(tc.sites, pat, tc.edge);
            c = s.type == "cuff";
            hw = s.sideXY_m(find(c, 1)) / 2;
            depth = tc.edge + R - sqrt(R^2 - hw^2);                 % edge indent + sagitta
            tc.verifyEqual(fit.centreIndent_m, depth * ones(nnz(c), 1), 'AbsTol', 1e-12);
            tc.verifyEqual(s.r_m(c), (R - depth + s.thickness_m(c)/2), 'AbsTol', 1e-12);
            tc.verifyEqual(fit.tiltFromRadial_deg, zeros(nnz(c), 1), 'AbsTol', 1e-6);
            % normal points into the nerve, as in sites.mat
            tc.verifyEqual(s.n(c, :), tc.sites.n(c, :), 'AbsTol', 1e-12);
        end

        function rotatedOvalFaceTangentAndSeated(tc)
            % on a rotated oval: face normal perpendicular to the outline tangent at the landmark, both side edges inside the nerve, the shallower one exactly edgeIndent deep
            pat = struct('nerveSemiAxes_m', [2.173e-3 1.347e-3], 'nerveRotation_deg', 23);
            [s, fit] = placeCuffContacts(tc.sites, pat, tc.edge);
            for k = find(s.type == "cuff")'
                ray = [cosd(s.theta_deg(k)) sind(s.theta_deg(k))];
                anchor = ray / sqrt(ellipseQ(ray, pat.nerveSemiAxes_m, pat.nerveRotation_deg));
                % tangent of the outline at the landmark, by central differences on the ellipse
                dth = 1e-6;
                rayP = ray * [cos(dth) sin(dth); -sin(dth) cos(dth)];
                rayM = ray * [cos(dth) -sin(dth); sin(dth) cos(dth)];
                pP = rayP / sqrt(ellipseQ(rayP, pat.nerveSemiAxes_m, pat.nerveRotation_deg));
                pM = rayM / sqrt(ellipseQ(rayM, pat.nerveSemiAxes_m, pat.nerveRotation_deg));
                tang = (pP - pM) / norm(pP - pM);
                nIn = s.n(k, 1:2);
                tc.verifyLessThan(abs(dot(nIn, tang)), 1e-6);       % face normal _|_ outline tangent
                tc.verifyLessThan(dot(nIn, anchor), 0);             % and pointing into the nerve
                % the two side edges of the inner face: edgeIndent below the native surface along n
                face = s.centre_m(k, 1:2) + nIn * s.thickness_m(k) / 2;
                for sgn = [-1 1]
                    e = face + sgn * s.sideXY_m(k) / 2 * s.w(k, 1:2);
                    tc.verifyLessThan(ellipseQ(e, pat.nerveSemiAxes_m, pat.nerveRotation_deg), 1);   % inside
                end
                j = fit.id == s.id(k);
                tc.verifyEqual(min(fit.edgeIndent_m(j, :)), tc.edge, 'AbsTol', 1e-12);
                % the shallower edge, moved back by edgeIndent along the outward normal, is on the
                % outline (Q = 1); the deeper one is still inside (Q < 1)
                qEdge = zeros(1, 2);
                for sgn = [-1 1]
                    e = face + sgn * s.sideXY_m(k) / 2 * s.w(k, 1:2) - tc.edge * nIn;
                    qEdge((sgn + 3) / 2) = ellipseQ(e, pat.nerveSemiAxes_m, pat.nerveRotation_deg);
                end
                tc.verifyEqual(max(qEdge), 1, 'AbsTol', 1e-9);
                % centre of the face deeper than the edges (convex outline)
                tc.verifyGreaterThan(fit.centreIndent_m(j), tc.edge);
            end
        end

        function shankContactsUntouched(tc)
            % placeCuffContacts must move only the cuff pads, never the intraneural (shank) contacts
            pat = struct('nerveSemiAxes_m', [2.173e-3 1.347e-3], 'nerveRotation_deg', 0);
            s = placeCuffContacts(tc.sites, pat, tc.edge);
            sh = tc.sites.type == "intraneural";
            tc.verifyEqual(s.centre_m(sh, :), tc.sites.centre_m(sh, :));
            tc.verifyEqual(s.n(sh, :), tc.sites.n(sh, :));
        end
    end
end
