function b = block_params()
%BLOCK_PARAMS Parameters of the conduction-block study (blocking the signal towards the amputated limb).
%   b = block_params()
%
%   Scenario (docs/design.md §5.3, tutor meeting 2026-11-23): the PROXIMAL site of the electrode
%   (shoulder side) stimulates a pattern of contacts to give a sensation; the action potentials (APs)
%   it launches towards the hand must not reach the reinnervated muscles (targeted muscle
%   reinnervation). The DISTAL site repeats the same pattern after a delay and/or with reversed
%   polarity (anodal block), so that the distal-going APs are cancelled. Choices: docs/decisions.md
%   D-12; walkthrough: docs/block.md.
%
%   The study re-uses the solved 3D model AS IT IS (results/main3d/<id>_3d.mat, the one shown in the
%   viewer): nothing is re-meshed or re-solved. The fibre model is MRG in NEURON (Python), driven by
%   the FEM potential sampled along straight fibres.
%
%   Used by
%     scripts/export_block_profiles.m   (MATLAB) samples the FEM potential along the fibres
%     scripts/run_block_study.py        (Python/NEURON) reads the same values from the JSON copy
%                                       written next to the profiles (results/block/<id>_block_params.json)
%
%   Tags: [source] literature / tutor; [group] group choice; [assumed] placeholder, to confirm or to
%   vary in the sensitivity analysis.
%   Units: SI (m, s, A, V). Frame: z = nerve axis, z = 0 at the centre of the device.

mm = 1e-3;  um = 1e-6;
b = struct();
b.set = "BLOCK-1";
b.patientId = "UA01";                       % the solved model in results/main3d/ (viewer data)

%% Sites ------------------------------------------------------------------------------------------
% The 48 contacts form two mirror-image halves (sections -4..-1 at z < 0, +1..+4 at z > 0, the same
% angles at +k and -k): these are the two sites. Which end faces the shoulder is not known; with the
% mirror symmetry it changes nothing, so z < 0 is called proximal.
b.proximalSide = "negative z";              % [assumed] D-12
% Distal contacts that are not solved are obtained from their solved mirror partner,
% V_k'(x, y, z) = V_k(x, y, -z). Exact for our model: the section is the same along z (extrusion,
% D-03) and the halves of the layout are mirror images (only the mesh differs: noise).
b.useMirror = true;                         % [group 2026-11-25] D-12

%% Fibres -----------------------------------------------------------------------------------------
% Straight fibres parallel to z (the model is an extrusion) on the lines of the viewer grid
% (scripts/export_viewer_data.m, 50 um), one every fibreStepXY_m, inside the endoneurium.
b.gridStepXY_m     = 50 * um;               % the viewer grid (must stay equal to export_viewer_data)
b.fibreStepXY_m    = 100 * um;              % [assumed] one fibre line every 2 viewer lines: ~500 fibres
                                            % in 10 fascicles (Zelechowski 2020: 240 fibres/mm^2 = 65 um)
b.fibreInset_m     = 10 * um;               % [assumed] keep lines at least this far inside the fascicle
                                            % circle (avoids lines on the perineurium itself)
b.profileStepZ_m   = 50 * um;               % sampling step of V along each fibre (MRG compartments are
                                            % 1 ... 200 um long; V is interpolated linearly between samples)
b.profileHalfSpan_m = 33 * mm;              % sample |z| <= 33 mm; the end planes (35 mm) are grounded

%% Fibre model: MRG (McIntyre, Richardson & Grill 2002, ModelDB 3810) ------------------------------
b.mrg.fibreDiam_um   = 10.0;                % [assumed] first runs: the MRG reference fibre (Aalpha/Abeta
                                            % range). Allowed: 5.7 7.3 8.7 10 11.5 12.8 14 15 16
b.mrg.fibreDiamSweep_um = [5.7 8.7 10.0 12.8 16.0];   % [assumed] diameters for the sweep
b.mrg.fibreLength_m  = 60 * mm;             % [assumed] centred on z = 0: ends 20 mm beyond the outer
                                            % stations, far from the stimulus (no end initiation)
b.mrg.nodeAlignment  = "proximalStation";   % [assumed] a node of Ranvier exactly at the z of the
                                            % proximal stimulating station (best case for excitation);
                                            % "random" = random offset per fibre (seeded)
b.mrg.seed           = 1;                   % rng seed for "random" alignment
b.mrg.temperature_C  = 37;                  % [source] MRG 2002
b.mrg.vInit_V        = -80e-3;              % [source] MRG 2002 resting potential
b.mrg.dt_s           = 5 * um;              % [source] MRG 2002 time step (5 us)
b.mrg.apThreshold_V  = 0;                   % an AP is an upward crossing of 0 mV at a node
b.mrg.detectInset_m  = 5 * mm;              % APs are counted at the nodes this far from each fibre end

%% Pulses -----------------------------------------------------------------------------------------
% Proximal (stimulation): monopolar cathodic rectangular pulse, as in Zelechowski 2020.
b.pulse.stimWidth_s   = 50 * um;            % [source] Zelechowski 2020 (candidate, see parameters.md)
b.pulse.stimOnset_s   = 0.5e-3;             % time of the proximal pulse after the start (rest settles)
% Distal (block): the same contacts mirrored; delay counted from the proximal onset.
b.pulse.blockWidth_s  = 50 * um;            % [assumed] collision: the same pulse as the proximal one
b.pulse.delaySweep_s  = (0:0.05:0.8) * 1e-3;   % [assumed] delays to test (arrival of a 10 um AP over
                                               % 20.6 mm is ~0.35 ms)
b.pulse.anodalWidth_s = 1.0e-3;             % [assumed] anodal block: long enough to cover the arrival
b.pulse.tAfter_s      = 3.0e-3;             % simulated time after the last pulse ends (propagation to
                                            % the fibre ends; slowest fibre 5.7 um ~ 30 m/s)

%% Patterns ---------------------------------------------------------------------------------------
% A pattern = proximal contacts and their weights (fraction of the pattern current, cathodic < 0
% for stimulation) and the fascicle it is meant to excite (the supervisor: "a configuration selected
% to excite a particular fascicle"). The distal pattern is the mirror partner of each contact
% (section -k -> +k). First runs: single contacts of station -4, the only one solved (contacts 1-6);
% targets = the fascicle each contact lies in (out.results.location: AS1 F2, AS2 F7, AS4 F1).
b.patterns = struct('name', {"AS1", "AS2", "AS4"}, ...
                    'contacts', {1, 2, 4}, ...
                    'weights',  {-1, -1, -1}, ...
                    'target',   {2, 7, 1});
% Speed: exact thresholds (bisection) only for the fibres of the target fascicle; every other fibre is
% simulated once at the stimulation current ("does it fire?"). --all-thresholds restores the full set.

%% Cathodal block on the distal site (supervisor, 2026-11-25) -------------------------------------
% "the same electrodes, on the other side, work as cathodes to cancel the signal". At the stimulation
% current this is the collision scenario (the distal AP still reaches the muscle, docs/block.md §6).
% A STRONG cathode can also block: the hyperpolarised regions on both sides of it (virtual anodes)
% stop the arriving AP and trap the cathode's own AP (seen at the proximal site: 400 uA, no AP out).
% Swept: distal cathodic amplitude (x I_stim) x onset delay x pulse width, on the target fascicle;
% the best settings are then re-checked on ALL fibres (a strong cathode excites other fascicles).
b.cathodal.factors  = [2 4 8 16 24 32 48 64];   % [assumed] distal amplitude / I_stim
b.cathodal.delays_s = (0:0.1:0.6) * 1e-3;       % [assumed] distal onset after the proximal onset
b.cathodal.widths_s = [50 200 500] * um;        % [assumed] distal pulse width
b.cathodal.nCheckAll = 3;                       % best settings re-checked on all fibres
b.searchStart_A = 5e-6;                     % threshold search: start, doubled until the fibre fires...
b.searchMax_A  = 2e-3;                      % ...up to this [A] (NaN threshold above it)
b.searchTol    = 0.02;                      % relative tolerance of the threshold bisection
b.recruitFraction = 0.10;                   % [source] Zelechowski 2020: threshold = 10% of the fibres
b.stimFactor   = 1.2;                       % [assumed] proximal amplitude = 1.2 x the 10% threshold of
                                            % the target fascicle (a sensation just above threshold)
b.anodalFactors = [1 2 4 8 16 32 64];       % [assumed] distal anodal amplitude / proximal amplitude
                                            % (each also run WITHOUT the proximal pulse: does the anode
                                            % itself excite the fibre through its virtual cathodes?)
b.nWorkers     = 0;                         % parallel NEURON processes (0 = number of cores - 1)
end
