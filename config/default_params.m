function params = default_params()
%DEFAULT_PARAMS All model parameters, in SI units, with their source.
%   params = default_params()
%
%   The single place where numbers are set (docs/architecture.md, principle 4).
%   Every value is also listed in docs/parameters.md with its source and status.
%   NaN = not chosen yet: the decision code (D-a ... D-h, see
%   docs/meetings/2026-11-10_kickoff.md and docs/decisions.md) says who settles it.
%   Code that needs a NaN value must stop with an error, never guess.
%
%   Owner: WP5 (structure); each WP owns the values of its own section.

mm = 1e-3;  um = 1e-6;
root = fileparts(fileparts(mfilename('fullpath')));      % project root

%% Geometry (WP1) -------------------------------------------------------------
params.geometry.level = "proximal forearm";       % D-h (proposal)
% Patient P00 = the tutor's fixed section (PW_lesson_5_11_students.m, section 1)
params.geometry.tutorSection.nerveRadius_m = 1.5 * mm;   % area 7.13 mm^2, Delgado-Martinez 2016 p. 6
params.geometry.tutorSection.tPeri_m       = 0.04 * mm;  % tutor's value (5.5-8 % of these diameters)
params.geometry.tutorSection.fascicles_m   = [ 0.00   0.00   0.36     % [x y rEndo] in mm -> m
                                               0.86   0.40   0.30
                                              -0.84   0.42   0.30
                                               0.10  -0.92   0.32
                                              -0.74  -0.58   0.24
                                               0.86  -0.52   0.24
                                               0.02   0.96   0.26 ] * mm;
% Random generator (sprint 3): statistics from Delgado-Martinez 2016, to be filled by WP1
params.geometry.generator = struct();

%% Domain and boundary conditions (WP3) ---------------------------------------
params.domain.outerRadiusFactor = 3;          % outer radius = 3 R (tutor's script); D-d
params.domain.length_m          = 40 * mm;    % proposal > device length (~21 mm); D-d
params.domain.outerBC           = "ground";   % "ground" | "insulated"; D-d

%% Materials (WP5) -------------------------------------------------------------
% Baseline = Davis et al. 2023, Table 1 (ASCENT defaults), as in the tutor's script; D-c.
% Tissue codes: 1 endoneurium, 2 perineurium, 3 epineurium, 4 surrounding medium (muscle)
params.materials.tissueName = ["endoneurium"; "perineurium"; "epineurium"; "medium"];
params.materials.sigmaT_Spm = [1/6;    1/1149; 1/6.3; 0.086];  % across the fibres [S/m]
params.materials.sigmaL_Spm = [1/1.75; 1/1149; 1/6.3; 0.35];   % along z [S/m]
params.materials.shankSigma_Spm = 6.67e-14;   % polyimide, Zelechowski 2020 (only if shank modelled)
params.materials.perineuriumMode = "explicit";% "explicit" | "scaledRhoT" | ...; D-f

%% Electrode (WP2) -------------------------------------------------------------
params.electrode.sitesFile   = fullfile(root, 'data', 'electrode', 'sites.mat');
params.electrode.siteModel   = "void";        % "void" (addVoid + current BC) | "volume"; D-e
params.electrode.activeFaces = "";            % "both" | "one" — ask the tutor (question 3)
params.electrode.withShank   = false;         % D-e, shank size unknown (question 2)

%% Mesh (WP2) ------------------------------------------------------------------
% Target element sizes per tissue: tutor's 2D values, to revise in the convergence study
params.mesh.h_m   = [50; 20; 120; 600] * um;  % endo, peri, epi, medium
params.mesh.hSite_m = NaN;                    % element size at the AS (WP2)
params.mesh.hGrad = 1.3;                      % growth rate (tutor's script)

%% Stimulation (WP4/WP5) -------------------------------------------------------
params.stim.unitCurrent_A = 1;                % each AS solved with 1 A, then scaled (linearity)
params.stim.pulseWidth_s  = 50 * um;         % Zelechowski 2020 (candidate)
params.stim.chargeRange_C = [0.5 60] * 1e-9;  % range explored in Zelechowski 2020 (candidate)

%% Activation (WP4) ------------------------------------------------------------
params.activation.internodeFactor     = 100;  % internode = 100 * D (MRG, as in Zelechowski 2020)
params.activation.recruitmentFraction = 0.10; % threshold = 10 % of the fascicle (Zelechowski 2020)
params.activation.fibreDiameters_m    = NaN;  % median-nerve distribution: to find (WP4)
params.activation.afThreshold         = NaN;  % activating-function threshold: to find (WP4)

%% Safety (WP5) ----------------------------------------------------------------
params.safety.shannonK      = NaN;            % Shannon criterion k: to find
params.safety.cicGold_Cpm2  = NaN;            % charge-injection capacity of gold: to find

%% Virtual patients (WP1) ------------------------------------------------------
params.patients.n     = NaN;                  % number of virtual patients: open decision
params.patients.seeds = [];
end
