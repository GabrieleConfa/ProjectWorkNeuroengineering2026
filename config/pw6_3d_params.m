function p = pw6_3d_params()
%PW6_3D_PARAMS Parameters of the 3D median-nerve model (upper arm, generated 10-fascicle section).
%   p = pw6_3d_params()
%
%   Starts from pw6_2d_params() (materials, contact size, perineurium thickness...) and adds or
%   overrides what the 3D model needs. EVERY modelling choice of the 3D pipeline is a field here,
%   so it can be changed without touching the code (src/main/main_pipeline_3d.m).
%   Choices and their reasons: docs/decisions.md (D-02 ... D-09); numbers and sources:
%   docs/parameters.md.
%
%   Tags used below:
%     [source]    value from the literature / PW6 group document / tutor
%     [group]     choice of the group (2026-11-13)
%     [assumed]   NOT from a source: a placeholder or a numerical choice, to confirm or to test
%                 in the sensitivity analysis
%
%   Units: SI (m, S/m, A). Frame: z = nerve axis, z = 0 at the centre of the device.

mm = 1e-3;  um = 1e-6;
p = pw6_2d_params();
p.set = "PW6-3D";

%% Virtual patient -------------------------------------------------------------------------------
p.patient.id   = "UA01";                % upper-arm patient no. 1 (results, meshes and caches are named after it)
p.patient.seed = 1;                     % rng seed of the fascicle layout (NaN = tutor's P00 section)
% Where the anatomy comes from (Gmsh mesher, src/mesh/gmshAnatomy.m; D-14). Running another patient
% or a toy experiment = changing id + seed (or source / file) here; nothing else in the code.
p.patient.source = "generator";         % "generator": our upper-arm generator makePatient (D-09)
                                        % "virtual":   Marco's virtual patient V<seed> (repulsion packing)
                                        % "file":      a JSON patient, p.patient.file (docs/patient_data.md)
p.patient.fascicleShape = "irregular";  % [tutor: fascicles are neither circles nor ellipses] "irregular":
                                        % random smooth shapes with the generator's centres and areas
                                        % (roundness 0.68 +/- 0.08, Pusnik 2025); "circle": the generator's circles
p.patient.file = "";                    % JSON file when source = "file"

%% Cross-section generator (upper arm) — src/geometry/generateFascicleLayout.m ----------------------
g = struct();
g.level               = "mid-upper arm";        % [group] PW6, Delgado-Martinez 2016
g.nerveArea_m2        = 9.20 * mm^2;            % [source] Pusnik 2025, histology mean (PW6) -> R = 1.711 mm
% Shape (D-10): ELLIPSE with the same area. roundness = minor / major axis (for an ellipse this is
% also the ImageJ "roundness" 4A/(pi*major^2)), so the semi-axes are R/sqrt(k) and R*sqrt(k):
% 2.173 x 1.347 mm for k = 0.62. The rigid cuff is shaped on this native oval (tight, D-05); the
% extrusion keeps the same section along the 7 cm (D-03). "circle" gives back the old round nerve.
g.nerveShape          = "ellipse";              % [group 2026-11-22] "ellipse" | "circle"
g.nerveRoundness      = 0.62;                   % [source, TO VERIFY] Pusnik 2025 Table 1, histology mean at
                                                % mid-arm, as cited on Marco's branch; the paper is not in
                                                % Papers/ yet
g.nerveRotation_deg   = 0;                      % [assumed] major axis along x of the electrode frame:
                                                % orientation of the device on the nerve not known (tutor)
g.nFascicles          = 10;                     % [source] Pusnik 2025 (PW6)
g.totalFascicleArea_m2 = 5.01 * mm^2;           % [source] Pusnik 2025: sum of fascicle (endoneurium) areas
g.largestFascicleArea_m2 = 1.49 * mm^2;         % [source] Pusnik 2025 / PW6: one large fascicle
                                                % (Delgado-Martinez 2016: cluster C, upper arm only)
% Sizes of the other fascicles: equivalent diameters drawn from a log-normal distribution with the
% between-fascicle coefficient of variation of the upper-arm median nerve, then all scaled together
% so that the total area is exactly totalFascicleArea_m2.
g.diameterCV          = 0.5;                    % [source] Brill & Tyler 2017: upper arm 0.6 +/- 0.3 mm
g.minEndoRadius_m     = 0.12 * mm;              % [assumed] smallest fascicle allowed (mesh); Delgado
                                                % cluster A = 0.13 mm^2 -> r = 0.20 mm
% Position: fascicles spread homogeneously over the section regardless of size
% (Delgado-Martinez 2016, Fig. 8): uniform random start, then overlaps are relaxed.
% Number of "inner" fascicles (not on the convex hull of the fascicle centres):
% ~ 0.58 (N - 6.08) (Brill & Tyler 2017) = 2.3 for N = 10. Layouts outside +/- tolerance are redrawn.
g.useInnerFascicleRule = true;                  % [source] Brill & Tyler 2017
g.innerFascicleTolerance = 1;                   % [assumed]
% Shank contacts vs fascicles (D-09):
%   "pierce"    the rigid shank goes straight through the tissue (as a transversal intraneural
%               electrode is implanted): a contact may lie across a fascicle boundary, where it
%               pierces the perineurium. Only nearly tangent configurations are redrawn (mesh).
%   "pushAside" every contact entirely inside one fascicle or entirely in the epineurium. Since the
%               fascicles are the same along z, all 16 shank footprints constrain the same section
%               and the layout becomes electrode-shaped: for comparison only.
g.contactRule         = "pierce";               % [assumed] "pierce" | "pushAside"
g.contactClearance_m  = 20 * um;                % [assumed] numerical (mesh quality)
g.minGapFascicles_m   = 25 * um;                % [assumed] epineurium between two perineuria (mesh)
g.minGapNerve_m       = 40 * um;                % [assumed] epineurium between a perineurium and the
                                                % nerve surface; must exceed the cuff indentation
g.maxRelaxIterations  = 4000;                   % numerical
g.maxRestarts         = 400;                    % numerical: new random start if relaxation fails
p.geometry.generator = g;
p.geometry.layout    = "generated";             % "generated" | "tutorP00"
p.geometry.level     = g.level;

%% Perineurium in 3D (D-06) ----------------------------------------------------------------------------
% THIN LAYER: the perineurium is not meshed as a volume. The boundary of each fascicle (circle of
% radius rEndo) is a surface across which the potential jumps; the current crossing it is
% J = (V_in - V_out) / (rho*t), with the sheet resistance rho*t = t_real / sigma_peri.
% Why: a meshed 40 um layer along 70 mm gave 27.7 million tetrahedra (test 2026-11-13); Pelot 2019
% found the thin layer ("contact impedance") gives thresholds similar to a meshed perineurium, and
% ASCENT (Davis 2023) uses it.
p.perineurium.representation  = "thinLayer";           % the only one implemented in 3D
p.perineurium.realThickness_m = p.geometry.tPeri_m;    % [source] 14 um (PW6 / Horn 2023)
% sheet conductance [S/m^2] = sigma_peri / t_real, computed in the pipeline from p.materials

%% Domain (D-02, D-04) --------------------------------------------------------------------------------
p.domain.nerveLength_m     = 70 * mm;           % [group]
p.domain.outerRadiusFactor = 3;                 % surrounding muscle up to 3 x the MAJOR semi-axis
                                                % (= 3 R for a circle, as in 2D)
p.domain.lateralBC         = "ground";          % outer surface of the muscle: "ground" (as in 2D)
p.domain.endBC             = "insulated";       % end planes z = +/- L/2: "ground" | "insulated" [group
                                                % 2026-11-27, D-17: grounded ends made every fibre fire
                                                % at its cut end under an anode (D-04 superseded)]

%% Cuff (D-05) ----------------------------------------------------------------------------------------
p.cuff.length_m    = 50 * mm;                   % [group] centred on the contacts (z = 0)
p.cuff.wallThickness_m = 0.5 * mm;              % [assumed] NOT KNOWN — ask the tutor. The wall is a
                                                % void (perfect insulator), so its thickness only sets
                                                % where the muscle starts outside the cuff.
p.cuff.fit = "tight";                           % [group] inner radius of the cuff = R, no saline gap

%% Contacts (D-07) ------------------------------------------------------------------------------------
% Contacts only (no shank body yet) [group]. Gold, 50 x 200 x 200 um (pw6_2d_params).
% Cuff contacts are FLAT and rigid [group] (D-11, src/electrode/placeCuffContacts.m): the face of
% each pad is tangent to the native oval at the pad's angle theta (sites.mat) and is pushed along that
% normal into the nerve until BOTH side edges of the face are cuffEdgeIndent_m below the native
% surface, so the whole 200 um face touches the tissue. The electrode does not bend; the nerve
% conforms: the tissue under the face is replaced by the rigid pad. The depth at the centre of the
% face is larger by the sagitta (w/2)^2/(2 rho), rho = local radius of curvature: 1.4 um on the flat
% side of the oval, 6.0 um at the ends of the major axis (2.9 um on the old circle).
p.electrode.cuffEdgeIndent_m = 7 * um;          % [assumed] gives the old 10 um at the centre on the
                                                % circle; a thinner sliver of tissue under the pad
                                                % would give bad elements
% Active contact: equipotential conductor carrying the current I. Modelled as a cell with a high
% conductivity and a uniform current source I / volume inside it (the source integrates to exactly I,
% the high conductivity makes the contact equipotential: the same physics as the 1 V + scaling of
% the 2D model, D-07).
p.electrode.contactSigma_Spm = 1e3;             % [assumed] numerical stand-in for gold (4.1e7 S/m
                                                % makes the matrix too ill-conditioned); a check in the
                                                % pipeline prints how equipotential the contact is
% Idle contacts: "insulating" (as in 2D) or "floating" (metal not connected = equipotential, no
% net current).
p.electrode.idleContacts  = "insulating";       % [group, as 2D]
p.electrode.insulatorSigma_Spm = 1e-7;          % [assumed] numerical stand-in for an insulator
p.electrode.pdmsSigma_Spm = 1e-7;               % [assumed] PDMS of the lead bodies and of the sleeve (Gmsh
                                                % model): unfilled PDMS is an insulator (~1e-12 S/m), 1e-7 is a
                                                % numerical stand-in. Marco's sensitivity: a filled, conductive
                                                % PDMS would lower the cuff-pad R by 12-76 % (ask the tutor).

%% Which contacts are solved (one 3D solve per contact) --------------------------------------------
p.solve3d.sections = [];                        % sections whose contacts are solved; [] = all 48.
                                                % All 48 since the Gmsh model (D-14): every station has its
                                                % own random lead insertion, so station +k is no longer the
                                                % mirror of -k (the mirror trick of D-12 does not apply).
                                                % (-4:-1 + mirror until 2026-11-26.)
p.solve3d.saveSingle = true;                    % store the potentials in single precision

%% 3D mesh (D-08, D-14) [assumed: first guess, to be fixed by the convergence study] ------------------
% Linear tetrahedra only: the system is assembled by src/fem/assembleConductionP1.m (needed for the
% perineurium thin layer, which the PDE Toolbox does not offer).
p.mesh3d.mesher    = "gmsh";                    % [group 2026-11-26, D-14] "gmsh": Marco's model (irregular
                                                % fascicles, bent PDMS leads, displaced fascicles, PDMS sleeve);
                                                % needs Python with the gmsh package (pip install gmsh; or set
                                                % NERVE_GMSH_PYTHON to a python.exe that has it).
                                                % "pdetoolbox": the previous model (straight shanks, circles).
p.mesh3d.cacheMode = "reuse";                   % "reuse": load the saved mesh of the same geometry from
                                                % results/mesh/ if there is one, else generate and save it;
                                                % "regenerate": always run the mesher (debugging the generator;
                                                % the previous file is kept as *_previous.mat). Gmsh never gives
                                                % the same mesh twice, so "reuse" is what keeps solved contacts valid.
p.mesh3d.fineHalfLength_m = 12 * mm;            % fine zone |z| < 12 mm (contacts span |z| <= 10.4 mm)
p.mesh3d.hContact_m = 30 * um;                  % on the faces of the contacts
p.mesh3d.hEndo_m    = 100 * um;                 % fascicle boundaries and endoneurium, fine zone
p.mesh3d.hNerve_m   = 150 * um;                 % nerve surface and epineurium, fine zone
p.mesh3d.hCoarse_m  = 500 * um;                 % all nerve faces outside the fine zone
p.mesh3d.hMax_m     = 1.2 * mm;                 % muscle, largest element
p.mesh3d.hMin_m     = 5 * um;
p.mesh3d.hGrad      = 1.4;
p.mesh3d.hEdge_m    = 5 * um;                   % [group 2026-11-26, D-14] element size on the contact edges
                                                % (0 = off). Marco: R converges with 5 um on the pad edges
                                                % (+2 % internal, +5 % cuff pads vs no edge refinement,
                                                % changes < 0.5 % below 5 um).
p.solve3d.solver    = "pcg";                    % [group 2026-11-26, D-14] "pcg" (ichol preconditioner; much
                                                % less memory; falls back to "direct" if it does not converge)
                                                % | "direct" (backslash, ~3-7 min and most of the 32 GB per contact)
p.solve3d.pcgTol    = 1e-9;
p.solve3d.assemblyBlockSize = 150000;           % numerical: tetrahedra per block when assembling the
                                                % stiffness matrix (bounds the temporary memory; from
                                                % Marco's branch, 2026-11-22)

%% Plots ----------------------------------------------------------------------------------------------
p.plot3d.afSpacing_m = 1 * mm;                  % node spacing of the activating-function preview:
                                                % internode = 100 D (MRG) for a 10 um fibre. Preview
                                                % only: phase 5 uses the fibre-diameter distribution.

%% Materials (3D) ------------------------------------------------------------------------------------
% sigmaT_Spm / sigmaL_Spm from pw6_2d_params: [endoneurium; perineurium; epineurium; muscle].
% In 3D the conductivity is the diagonal tensor diag(sigmaT, sigmaT, sigmaL) (z = fibre direction).
p.materials.tissueName3d = ["endoneurium"; "perineurium"; "epineurium"; "muscle"; "contact"];
end
