function sites = loadActiveSites(filePath, nerveRadius_m)
%LOADACTIVESITES Load the 48 active sites (AS) of the provided electrode, in SI units.
%   sites = loadActiveSites(filePath, nerveRadius_m)
%
%   Reads the matrix "sites" (48 x 10, mm and degrees) provided by the tutor
%   (see data/electrode/README.md for its origin and column meaning), converts
%   it to SI units and adds the geometric vectors needed to build each site as
%   a box in the 3D model.
%
%   Inputs
%     filePath       path of the .mat file with the variable "sites"
%                    (normally data/electrode/sites.mat)
%     nerveRadius_m  nerve radius R [m]. Needed only for the cuff sites, whose
%                    centre is at r = R + thickness/2 so that their face rests
%                    on the nerve surface (instruction in load_active_sites.m).
%
%   Output: table, one row per site, columns
%     id            electrode number (1...48)
%     type          "intraneural" or "cuff" (categorical)
%     section       section number, -4...-1 (z<0) and +1...+4 (z>0)
%     z_m           position of the centre along the nerve axis [m]
%     r_m           distance of the centre from the nerve axis [m]
%     theta_deg     angular position of the centre, from x, counterclockwise [deg]
%     phi_deg       direction of the normal to the square face [deg]
%     thickness_m   thickness of the box, along the normal n [m]
%     sideZ_m       side of the face parallel to the nerve axis [m]
%     sideXY_m      side of the face in the x-y plane, along w [m]
%     centre_m      [x y z] of the box centre [m]           (1x3 per row)
%     n             unit normal, [cos(phi) sin(phi) 0]       (1x3 per row)
%     w             unit in-plane side, [-sin(phi) cos(phi) 0] (1x3 per row)
%
%   Reference frame (from the tutor): z = nerve axis, origin at the centre of
%   the device; x = r cos(theta), y = r sin(theta).

arguments
    filePath (1,:) char
    nerveRadius_m (1,1) double {mustBePositive}
end

mm = 1e-3;                                   % mm -> m
S  = load(filePath, "sites");
M  = S.sites;
assert(size(M, 2) == 10, 'loadActiveSites:format', ...
       'Expected a matrix with 10 columns, got %d.', size(M, 2));

isCuff      = M(:,2) == 2;
thickness_m = M(:,8) * mm;

r_m          = M(:,5) * mm;
r_m(isCuff)  = nerveRadius_m + thickness_m(isCuff) / 2;   % cuff face on the nerve surface
assert(~any(isnan(r_m)), 'loadActiveSites:nanRadius', ...
       'Some intraneural sites have no radius in the data file.');

typeName          = repmat("intraneural", size(M,1), 1);
typeName(isCuff)  = "cuff";

theta = M(:,6);   phi = M(:,7);
centre_m = [r_m .* cosd(theta), r_m .* sind(theta), M(:,4) * mm];
n        = [cosd(phi),  sind(phi), zeros(size(phi))];
w        = [-sind(phi), cosd(phi), zeros(size(phi))];

sites = table(M(:,1), categorical(typeName), M(:,3), M(:,4) * mm, r_m, theta, phi, ...
              thickness_m, M(:,9) * mm, M(:,10) * mm, centre_m, n, w, ...
              'VariableNames', {'id', 'type', 'section', 'z_m', 'r_m', 'theta_deg', ...
                                'phi_deg', 'thickness_m', 'sideZ_m', 'sideXY_m', ...
                                'centre_m', 'n', 'w'});
end
