function T = functionCatalog() % Functional groups carried by the median nerve at mid-arm (what a fibre controls or reports).
%% Purpose
% A fascicle at mid-arm is MIXED (Delgado-Martinez 2016): it carries fibres of several groups. The model therefore stores, for
% every fascicle, a composition: the fraction of its fibres that belongs to each group below. A patient-specific measurement can
% set the composition (for example one group = 1); a virtual patient draws it at random (assignFunctions).
% Columns: id (short key), name, kind (motor or sensory), target (muscles or skin), exit_pctForearm (where the branch leaves
% the trunk, % of forearm length, 0 = wrist, 100 = elbow; NaN = not in the forearm or not given), priorWeight (share of the motor
% or sensory fibres of the nerve: ASSUMED, uniform, no data), colour (RGB, figures and viewer), evidence.
% No inputs. Output T: table with one row per group (7 rows); its row order is the order of every "composition" vector.
id =["PT";"FCRPL";"FDS";"AIN";"THENAR";"PCB";"DIGITAL"]; % Keys used in patient files.
name = ["Pronator teres";"Flexor carpi radialis and palmaris longus";"Flexor digitorum superficialis";"Anterior interosseous (FDP, FPL, pronator quadratus)"; ...
        "Thenar and lumbrical muscles (recurrent branch)";"Palmar cutaneous branch";"Digital sensory (thumb, index, middle, half ring)"];
kind = ["motor";"motor";"motor";"motor";"motor";"sensory";"sensory"];
target = ["pronation";"wrist flexion";"finger flexion (PIP)";"distal finger and thumb flexion, pronation";"thumb opposition, finger dexterity";"palm skin";"finger skin"];
exit_pctForearm = [98.4;NaN;45.8;81.0;NaN;16.5;NaN]; % Delgado-Martinez 2016: PT upper branch 98.4, interosseous 81.0, FDS 45.8, palmar cutaneous 16.5 (n = 3 to 5 nerves). FCR/PL: P0.5 at about 84 and 110 percent (not an exit point).
priorWeight = [0.2;0.2;0.2;0.2;0.2;0.2;0.8]; % ASSUMED. Motor groups: uniform. Sensory: the digital fibres are taken as the larger part (no data). Normalised within each kind.
colour = [0.85 0.37 0.01;0.80 0.47 0.74;0.00 0.45 0.70;0.94 0.89 0.26;0.00 0.62 0.45;0.34 0.71 0.91;0.60 0.60 0.60];
evidence = ["Delgado-Martinez 2016 (exit points)";"Delgado-Martinez 2016 (branching probability)";"Delgado-Martinez 2016";"Delgado-Martinez 2016; MR neurography 2014 (dorsal and lateral in the upper arm)"; ...
            "general anatomy (not in the sources read)";"Delgado-Martinez 2016";"general anatomy (not in the sources read)"];
T = table(id,name,kind,target,exit_pctForearm,priorWeight,colour,evidence); % One row per group.
end
