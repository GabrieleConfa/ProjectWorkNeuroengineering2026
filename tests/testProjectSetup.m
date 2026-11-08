classdef testProjectSetup < matlab.unittest.TestCase
%TESTPROJECTSETUP Checks that the shared configuration loads and is in SI units.
%   Run with:  setup_paths; runtests('tests/testProjectSetup.m')

    methods (Test)
        function paramsLoad(tc)
            p = default_params();
            for f = ["geometry" "domain" "materials" "electrode" "mesh" ...
                     "stim" "activation" "safety" "patients"]
                tc.verifyTrue(isfield(p, f), "missing params." + f);
            end
        end

        function sitesFileExists(tc)
            p = default_params();
            tc.verifyTrue(isfile(p.electrode.sitesFile), p.electrode.sitesFile);
        end

        function lengthsAreInMetres(tc)
            % catches a value typed in mm by mistake (nerve radius ~1e-3 m, not ~1)
            p = default_params();
            tc.verifyLessThan(p.geometry.tutorSection.nerveRadius_m, 0.01);
            tc.verifyLessThan(max(p.geometry.tutorSection.fascicles_m(:)), 0.01);
            tc.verifyLessThan(p.domain.length_m, 0.2);
            tc.verifyLessThan(max(p.mesh.h_m), 0.01);
        end

        function tutorFasciclesInsideNerve(tc)
            % every fascicle (endoneurium + perineurium) lies inside the nerve
            s = default_params().geometry.tutorSection;
            F = s.fascicles_m;
            outer = hypot(F(:,1), F(:,2)) + F(:,3) + s.tPeri_m;
            tc.verifyLessThan(outer, s.nerveRadius_m * ones(size(outer)));
        end

        function anisotropyAlongZ(tc)
            % endoneurium conducts better along the fibres than across them
            m = default_params().materials;
            tc.verifyGreaterThan(m.sigmaL_Spm(1), m.sigmaT_Spm(1));
        end
    end
end
