function s = insertionParams() % Parameters of the lead-insertion Monte Carlo (SI units). Every value is tagged.
% Tags: [GIVEN] told by the tutor, [LIT] from a paper (order of magnitude, see docs/algorithms_and_plan.md),
%       [ASSUMED] our choice, NOT calibrated: results are indicative, not probabilities.
%% 1. Lead and tip
s.tipDiameter_m      = 100e-6;   % [GIVEN] tip of 100 um, interpreted as a diameter.
s.leadYoung_Pa       = 2e6;      % [LIT] PDMS (Sylgard 184, 10:1) tensile modulus 1.3-3 MPa (arXiv:2404.16960). Reported only: motion is displacement-controlled.
s.stepLength_m       = 10e-6;    % [ASSUMED] numerical step of the advancing tip.
%% 2. How the clinician inserts it
s.entryAngleRange_deg = 0;       % [ASSUMED] uniform +/- range of the push direction around the nominal axis. 0 since
                                 % 2026-11-26 (tutor: the exit is on the line of the entry); Marco's force model used 10.
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
%% 4b. Pierce model (added on main, 2026-11-26)
s.pierceModel     = "probability"; % [group 2026-11-26, D-14] production model of the Gmsh mesh (insertionScenario).
                                 % "force": Marco's model above (forces compared at each contact: slide / displace / pierce).
                                 % "probability": the tutor's rule (meeting 2026-11-23) with the angle added (user
                                 % 2026-11-26): when the tip touches a fascicle it pierces it with probability
                                 % P = 1/(1 + (r50/r)^n) * cos(gamma)^m (src/insertion/pierceProbability.m), r = sqrt(A/pi),
                                 % gamma = angle to the perineurium normal at the first contact; otherwise it bends around
                                 % it (no displacement). Use it with returnToLine = true and no random entry angle/offset,
                                 % so that the exit stays on the line of the entry.
s.pierceR50_m     = 0.38e-3;     % [ASSUMED] radius with a 50 % head-on probability. With n = 6, on UA01 (radii 0.17-0.69 mm):
                                 % F10 0.17 mm -> 1 %, F7 0.28 mm -> 13 %, F3 0.43 mm -> 68 %, F1 0.69 mm -> 97 % (head-on):
                                 % the tutor's "small fascicles very low, large ones high". To ask the tutor.
s.pierceSteepness = 6;           % [ASSUMED] exponent n: how sharp the step from small to large fascicles is.
s.pierceAngleExponent = 2;       % [ASSUMED] exponent m of cos(gamma)^m: P halves at 45 deg, quarters at 60 deg.
s.returnToLine    = true;        % steer the tip back to the line of the entry in free motion (tutor's rule; Marco: false).
s.returnLookahead_m = 0.3e-3;    % [ASSUMED] look-ahead distance of that steering (pure pursuit): smaller = sharper return.
%% 5. Monte Carlo
s.nTrials = 100;                 % Insertions simulated per station.
s.seed    = 20261124;            % Reproducible random streams.
end
