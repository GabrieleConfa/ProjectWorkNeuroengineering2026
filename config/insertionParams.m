function s = insertionParams() % Parameters of the lead-insertion Monte Carlo (SI units). Every value is tagged.
% Tags: [GIVEN] told by the tutor, [LIT] from a paper (order of magnitude, see docs/algorithms_and_plan.md),
%       [ASSUMED] our choice, NOT calibrated: results are indicative, not probabilities.
%% 1. Lead and tip
s.tipDiameter_m      = 100e-6;   % [GIVEN] tip of 100 um, interpreted as a diameter.
s.leadYoung_Pa       = 2e6;      % [LIT] PDMS (Sylgard 184, 10:1) tensile modulus 1.3-3 MPa (arXiv:2404.16960). Reported only: motion is displacement-controlled.
s.stepLength_m       = 10e-6;    % [ASSUMED] numerical step of the advancing tip.
%% 2. How the clinician inserts it
s.entryAngleRange_deg = 10;      % [ASSUMED] uniform +/- range of the push direction around the nominal axis.
s.entryOffsetSD_m     = 0.15e-3; % [ASSUMED] sideways scatter of the entry point (standard deviation).
s.pushForceMedian_N   = 20e-3;   % [ASSUMED] normal force a hand can deliver at the tip (log-normal median).
s.pushForceSigma      = 0.5;     % [ASSUMED] log-normal shape parameter of the hand force.
%% 3. Resistances compared at every encounter
s.punctureForceMedian_N = 2e-3;  % [LIT] > 2 mN to puncture rabbit nerve with a tungsten needle (Jensen 2001). Transferred to a 100 um PDMS tip: [ASSUMED].
s.punctureForceSigma    = 0.4;   % [ASSUMED] log-normal variability of the puncture force.
s.fascicleStiffness_Npm = 20;    % [ASSUMED] force per metre of fascicle displacement in the epineurium (F_displace = k * displacement).
s.frictionCoeff         = 0.5;   % [ASSUMED] the tip slides along a perineurium when tan(angle to the normal) exceeds this value.
s.fascicleGap_m         = 25e-6; % [ASSUMED] minimum gap kept between displaced fascicles (as in the layout generator).
%% 4. Shape statistics of the fascicles
s.roundnessMean = 0.68;          % [LIT] fascicle roundness 4A/(pi*major^2), Pusnik 2025, histology.
s.roundnessSD   = 0.08;          % [LIT] same source.
s.harmonicAmplitude = 0.07;      % [ASSUMED] amplitude of the irregularity (radial Fourier harmonics 3 to 5).
%% 5. Monte Carlo
s.nTrials = 100;                 % Insertions simulated per station.
s.seed    = 20261124;            % Reproducible random streams.
end
