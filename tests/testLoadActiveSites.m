classdef testLoadActiveSites < matlab.unittest.TestCase
%TESTLOADACTIVESITES Checks on the provided active-site layout and its loader.
%   Run with:  setup_paths; runtests('tests/testLoadActiveSites.m')
%   These checks encode what load_active_sites.m (tutor) states about the
%   device, so a corrupted or different data file is caught immediately.

    properties
        sites
        R_m = 1.5e-3     % test nerve radius [m] (value of the tutor's 2D script)
    end

    methods (TestClassSetup)
        function loadData(tc)
            root = fileparts(fileparts(mfilename('fullpath')));
            tc.sites = loadActiveSites(fullfile(root, 'data', 'electrode', 'sites.mat'), tc.R_m);
        end
    end

    methods (Test)
        function countsPerSection(tc)
            % 48 sites, 8 sections, 4 intraneural + 2 cuff in each section
            s = tc.sites;
            tc.verifyEqual(height(s), 48);
            tc.verifyEqual(sort(unique(s.section))', [-4:-1, 1:4]);
            for k = [-4:-1, 1:4]
                in = s.section == k;
                tc.verifyEqual(nnz(in & s.type == "intraneural"), 4);
                tc.verifyEqual(nnz(in & s.type == "cuff"), 2);
            end
        end

        function sectionPositions(tc)
            % sections 2 mm apart, sections -1/+1 at 4.3 mm from z = 0
            zSec = unique(abs(tc.sites.z_m));
            tc.verifyEqual(zSec', [4.3 6.3 8.3 10.3] * 1e-3, 'AbsTol', 1e-12);
        end

        function mirrorSymmetry(tc)
            % section +k is the mirror image of section -k with respect to z = 0
            s = tc.sites;
            for k = 1:4
                a = sortrows([s.centre_m(s.section == -k, :), s.phi_deg(s.section == -k)]);
                b = sortrows([s.centre_m(s.section ==  k, :), s.phi_deg(s.section ==  k)]);
                b(:,3) = -b(:,3);
                tc.verifyEqual(sortrows(a), sortrows(b), 'AbsTol', 1e-12);
            end
        end

        function radii(tc)
            % intraneural sites at 0.3 and 0.9 mm; cuff sites at R + thickness/2
            s = tc.sites;
            intra = s.type == "intraneural";
            tc.verifyEqual(sort(unique(round(s.r_m(intra) * 1e6)))', [300 900]);
            cuff = s.type == "cuff";
            tc.verifyEqual(s.r_m(cuff), tc.R_m + s.thickness_m(cuff) / 2, 'AbsTol', 1e-12);
        end

        function orthonormalFrame(tc)
            % n and w are unit vectors, orthogonal to each other and to z
            s = tc.sites;
            tc.verifyEqual(vecnorm(s.n, 2, 2), ones(48,1), 'AbsTol', 1e-12);
            tc.verifyEqual(vecnorm(s.w, 2, 2), ones(48,1), 'AbsTol', 1e-12);
            tc.verifyEqual(sum(s.n .* s.w, 2), zeros(48,1), 'AbsTol', 1e-12);
            tc.verifyEqual(s.n(:,3), zeros(48,1));
        end

        function cuffNormalsPointInward(tc)
            % cuff sites: the normal points towards the nerve axis
            s = tc.sites(tc.sites.type == "cuff", :);
            radial = s.centre_m(:,1:2) ./ vecnorm(s.centre_m(:,1:2), 2, 2);
            tc.verifyEqual(sum(s.n(:,1:2) .* radial, 2), -ones(height(s),1), 'AbsTol', 1e-12);
        end

        function intraneuralFacesAlongShank(tc)
            % intraneural sites: the face lies along the shank, i.e. w is radial
            s = tc.sites(tc.sites.type == "intraneural", :);
            radial = s.centre_m(:,1:2) ./ vecnorm(s.centre_m(:,1:2), 2, 2);
            tc.verifyEqual(abs(sum(s.w(:,1:2) .* radial, 2)), ones(height(s),1), 'AbsTol', 1e-12);
        end
    end
end
