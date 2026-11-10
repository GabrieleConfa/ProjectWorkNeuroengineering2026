function p = pw6_2d_params()
%PW6_2D_PARAMS Parameters of the 2D median-nerve model selected by the group (PW6), in SI units.
%   p = pw6_2d_params()
%
%   Source of every selected value: report/PW6_model_choices_and_rationale.pdf (group document,
%   2026-11-10), which also gives the reason and the literature reference of each choice.
%   Values NOT given in that document are marked [NOT IN PW6] with what is used instead; they
%   must be confirmed by the group.
%
%   Starts from default_params() so that the fields not touched here (domain, stimulation,
%   activation, safety...) keep their project defaults.

mm = 1e-3;  um = 1e-6;
p = default_params();
p.set = "PW6";

%% Anatomy ------------------------------------------------------------------------------------
% Selected in PW6 (kept here for the record and for the next geometry version):
p.geometry.level = "mid-upper arm, before the forearm motor branches";   % Delgado-Martinez 2016
nerveArea_m2 = 9.20 * mm^2;                       % Pusnik 2025, histological cross-sections (HCS) mean
p.geometry.pw6.nerveRadius_m = sqrt(nerveArea_m2 / pi);   % equivalent circle, 1.711 mm
p.geometry.pw6.nFascicles     = 10;               % Pusnik 2025, HCS mean count
p.geometry.pw6.areaTotal_m2   = 5.01 * mm^2;      % total fascicle area
p.geometry.pw6.areaLargest_m2 = 1.49 * mm^2;      % one large fascicle, the other nine equal
% [NOT IN PW6] the fascicle POSITIONS. The PW6 sizes do not fit in the tutor's 7 positions,
% so for the current runs the group chose (2026-11-10) to keep the PREVIOUS geometry:
% the tutor's 7-fascicle section (patient P00, R = 1.5 mm) with the central fascicle enlarged
% to an endoneurium radius of 0.42 mm, so that the two inner contacts lie inside it.
p.geometry.layout = "tutorP00";
p.geometry.centralREndo_m = 0.42 * mm;

% Perineurium thickness: the PW6 mesh note states "the perineurium is 14 um thick"
% (consistent with Horn 2023, ~14 um measured on canine fascicles). Applied to the P00 layout.
p.geometry.tPeri_m = 14 * um;

%% Electrode ----------------------------------------------------------------------------------
p.electrode.zCut_m = -10.3 * mm;                  % device station drawn in 2D (section -4)
p.electrode.thicknessOverride_m = 50 * um;        % gold contacts 50 x 200 x 200 um (slides);
                                                  % overrides the 5 um of sites.mat
p.electrode.material = "gold";
p.electrode.cuffPlacement = "surface";            % cuff centre at R + thickness/2, no squeeze
p.electrode.siteModel = "ideal";                  % gold = ideal contact boundary (equipotential),
                                                  % no bulk gold conductivity in the mesh
% [NOT IN PW6] idle contacts: insulating holes. Stimulation: each contact monopolar against
% the grounded outer boundary, plus bipolar pairs inside the station.

%% Stimulation [NOT IN PW6] ------------------------------------------------------------------
p.stim.currentPerLength_Apm = 1;                  % 2D "unit current": 1 A per metre of nerve
                                                  % (linear problem: results scale with it)
% Bipolar pairs, one row per configuration: [cathode anode] as the POSITION of the contact
% inside the station (1-4 shank contacts in file order, 5-6 cuff contacts). Station -4:
% 1 = r 0.9 mm at +x, 2 = r 0.3 mm at +x, 3 = r 0.3 mm at -x, 4 = r 0.9 mm at -x,
% 5 = cuff at +y, 6 = cuff at -y.
p.stim.bipolarPairs = [2 3                        % both inside the central fascicle
                       2 1                        % inner (in a fascicle) vs outer, same side
                       5 6];                      % the two cuff contacts, across the nerve

%% Tissue conductivities [S/m], order: endoneurium, perineurium, epineurium, surrounding muscle
p.materials.tissueName = ["endoneurium"; "perineurium"; "epineurium"; "muscle"];
p.materials.sigmaT_Spm = [0.0826; 0.00088; 0.0826; 0.13];  % Zelechowski 2020 / Pelot 2019; muscle Kangasmaa 2025
p.materials.sigmaL_Spm = [0.571;  0.00088; 0.0826; 0.56];  % kept for 3D; not used in the 2D solve

%% Mesh [m] -----------------------------------------------------------------------------------
p.mesh.h_m(2)   = 5 * um;                         % hPeri = 0.005 mm (several elements across 14 um)
p.mesh.hSite_m  = 10 * um;                        % hElec = 0.01 mm on the contact edges
p.mesh.hMin_m   = 2.5 * um;                       % [NOT IN PW6] minimum size = hPeri / 2
% [NOT IN PW6] endoneurium / epineurium / medium sizes: project defaults (50, 120, 600 um).

%% Domain -------------------------------------------------------------------------------------
% [NOT IN PW6] outer radius of the medium and its boundary condition: project defaults
% (3 R, grounded).
end
