function p = quicklook3DParams() % Collect the editable inputs for the 3D geometry preview.
%% 1. Provisional nerve section
p = pw6_2d_params(); % Reuse the group settings and tutor's contact file.
p.domain.length_m = 60e-3; % Show the full 6 cm segment discussed with the tutor.
p.geometry.radius_m = p.geometry.pw6.nerveRadius_m; % Use the area-equivalent PW6 circle for now.
% These illustrative positions allow both a straight crossing and a clear route at every station.
fascicles_mm = [0.000,  0.000, 0.38;  % F1, central demonstration obstacle.
                1.245,  0.345, 0.21;  % F2, peripheral.
                0.495,  1.245, 0.21;  % F3, peripheral.
               -0.430,  1.245, 0.21;  % F4, peripheral.
               -1.155,  0.420, 0.21;  % F5, peripheral.
               -0.780,  0.020, 0.21;  % F6, peripheral.
               -1.255, -0.505, 0.21;  % F7, peripheral.
               -0.355, -1.305, 0.21;  % F8, peripheral.
                0.670, -1.005, 0.21;  % F9, peripheral.
                1.245, -0.480, 0.21]; % F10, peripheral.
p.geometry.fascicles_m = 1e-3 * fascicles_mm; % Convert [x y endoneurium radius] from mm to metres.
%% 2. Cuff and strip geometry
p.electrode.cuffIndentation_m = 0.12e-3; % Prescribe a flat contact depth for the first geometry check.
p.electrode.cuffTransition_m = 0.40e-3; % Blend back to the circular nerve beyond each cuff pad edge.
p.electrode.maxDeflection_m = 0.8e-3; % Maximum lateral shift explored by the geometric search.
%% 3. Provisional strip mechanics for a puncture feasibility screen
p.electrode.stripWidth_m = 0.12e-3; % Assume a flexible strip narrower than its 0.20 mm active pads.
p.electrode.stripThickness_m = 10e-6; % Assume a thin flexible substrate pending the device specification.
p.electrode.stripYoungModulus_Pa = 2.5e9; % Use published photosensitive-polyimide modulus as a material proxy.
p.electrode.unsupportedLength_m = 1e-3; % Assume one millimetre can bend near a fascicle contact.
p.electrode.perineuriumPunctureForce_N = 40e-3; % Use the measured rabbit perineurium probe force as a provisional comparison.
end % Return geometry and electrode settings in SI units.
