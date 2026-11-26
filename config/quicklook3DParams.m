function p = quicklook3DParams(opts) % Collect the editable inputs for the 3D geometry preview. opts (optional): patientFile, virtualSeed, functionScenario, functionSeed.
if nargin<1, opts = struct(); end % Defaults: virtual patient 1.
%% 1. Provisional nerve section
p = pw6_2d_params(); % Reuse the group settings and tutor's contact file.
p.domain.length_m = 60e-3; % Show the full 6 cm segment discussed with the tutor.
p.geometry.radius_m = p.geometry.pw6.nerveRadius_m; % Retain a scalar scale for the provisional circle.
p.geometry.roundness=0.62; % Pusnik 2025, Table 1, histological mean at mid-arm.
p.geometry.nerveSemiAxes_m=p.geometry.radius_m*[1/sqrt(p.geometry.roundness),sqrt(p.geometry.roundness)]; % Preserve the measured 9.20 mm2 section.
p.geometry.nerveRotation_deg=0; % Provisional device-to-nerve orientation.
p.geometry.layoutSeed=1; % Synthetic coordinates; individual fascicle positions were not measured.
p.geometry.layoutTemplate=[ ... % Synthetic reference: normalized centres and relative fascicle areas.
    -0.2280128549 -0.4214593498 1.4900000000; ... % F1, simulated position.
    -0.4434229617 0.3119792249 0.1276734544; ... % F2, simulated position.
    -0.7915181022 0.2004375260 0.3404542581; ... % F3, simulated position.
    -0.7522162956 -0.3760097121 0.1946951764; ... % F4, simulated position.
    0.5893623139 -0.2436962103 1.1899043082; ... % F5, simulated position.
    0.2164414934 -0.7932066252 0.1179122677; ... % F6, simulated position.
    0.6174268395 0.5438026978 0.1872820740; ... % F7, simulated position.
    -0.3162605864 0.7383269006 0.1700136156; ... % F8, simulated position.
    0.1375707486 0.5124654774 1.1170499756; ... % F9, simulated position.
    -0.6142282775 0.6225495706 0.0750148700]; % F10, simulated position.
p.geometry.fascicles_m=midArmLayout(p); % Match the 5.01 mm2 total and 1.49 mm2 largest fascicle (circular equivalents).
p.geometry.mode="insertion"; % "insertion": random irregular fascicles and random lead insertions; "legacy": circular fascicles and the approved routes.
%% 2. Cuff and strip geometry
p.mechanics.mu_Pa=12.9e3; % Main et al. (2011), Table 1, whole median nerve at the wrist.
p.mechanics.alpha=6.5; % First-order Ogden exponent with the published mu/alpha^2 convention.
p.electrode.cuffMeasuredDepth_m = zeros(0,2); % Optional rows [external site ID, measured centre indentation in metres].
p.electrode.cuffExtraIndentation_m = 3e-6; % Trial closure per pad beyond full seating (assumed, adjustable).
p.electrode.cuffTransition_m = 0.40e-3; % Blend along the nerve beyond each cuff pad edge.
p.electrode.cuffLateralTransition_m = 0.12e-3; % Blend around the 0.20 mm pad width in the cross-section.
p.electrode.maxDeflection_m = max(p.geometry.nerveSemiAxes_m); % Maximum lateral shift explored by the geometric search.
p.electrode.routingModel = "beam"; % Use reduced-order elastic equilibrium with perineurial contact.
p.electrode.foundation_Npm2 = 250; % Provisional lateral epineurial support per unit length per displacement.
p.electrode.reverseEntrySections = 1:4; % Provisional winding layout: insert from the opposite end of each repeated axis.
p.electrode.entryAngleRange_deg = 5; % Illustrative symmetric insertion tolerance in degrees; replace with measured data.
p.electrode.initialBow_m = 0.045e-3; % Provisional bound on the unloaded strip imperfection.
p.electrode.appliedInsertionForce_N = NaN; % Supply a measured insertion load to enable the conditional rupture scenario.
p.electrode.calibratedPunctureForce_N = NaN; % Supply a device-specific normal puncture force when measured.
p.electrode.routeSeed = 20261120; % Editable seed for provisional independent insertion trajectories.
p.electrode.tipOverhang_m = 0.05e-3; % Provisional distance beyond the distal pad edge along the insertion axis (legacy mode only).
%% 3. Local mesh controls
p.meshPreview.hMax_m = 0.18e-3; % Set the first tetrahedral target size for a single contact slab.
p.meshPreview.hMin_m = 0.025e-3; % Resolve the local contact without implying a convergence study.
%% 4. Provisional strip mechanics for a puncture feasibility screen
p.electrode.stripWidth_m = 0.12e-3; % Assume a flexible strip narrower than its 0.20 mm active pads.
p.electrode.stripThickness_m = 10e-6; % Assume a thin flexible substrate pending the device specification.
p.electrode.stripYoungModulus_Pa = 2.5e9; % Use published photosensitive-polyimide modulus as a material proxy.
p.electrode.unsupportedLength_m = 1e-3; % Assume one millimetre can bend near a fascicle contact.
p.electrode.perineuriumPunctureForce_N = 40e-3; % Use the measured rabbit perineurium probe force as a provisional comparison.
%% 5. Patient: virtual V001 (default) or a specific patient added as a JSON file (see docs/patient_data.md)
p.patient = struct('name',"",'file',"",'virtualSeed',1,'functionScenario',"none",'functionSeed',1); % Default: reference virtual patient V001 WITHOUT function assignment ("none"): the functions of the fascicles are not known. Hypotheses: "random", "forearmRule", "mrAIN". name = file in data/patients/personal; file = any JSON path.
for name = string(fieldnames(opts))', if name=="patientFile", p.patient.file = string(opts.patientFile); else, p.patient.(name) = opts.(name); end, end % Overrides.
p.geometry.patient = []; p.geometry.scenario = [];
if p.geometry.mode=="insertion"
    p.geometry.patient = buildPatient(p); % Section, fascicle outlines, functions and motor patches of the patient used in this run.
    if ~(p.geometry.patient.source=="virtual" && p.geometry.patient.seed==1) % Legacy circular representation (used by some older helpers) follows the patient.
        o = p.geometry.patient.fascicles; p.geometry.fascicles_m = cell2mat(arrayfun(@(f)[mean(f.outline,1),sqrt(f.area_m2/pi)*[1 1],0],o(:),'UniformOutput',false));
    end
    p.geometry.scenario = insertionScenario(p); % One random lead insertion per station on this patient, shared by figures, mesh and viewer.
end
end % Return geometry and electrode settings in SI units.
