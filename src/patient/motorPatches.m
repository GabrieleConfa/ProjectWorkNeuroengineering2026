function patches = motorPatches(outline, motorFraction, rs, margin_m) % Random "patches" of motor fibres inside one fascicle.
%% Purpose
% Delgado-Martinez 2016: motor axons (about one eighth of all axons) are concentrated in patches inside fascicles, and the
% position of the patches does not depend on the size or position of the fascicle. Here 1 to 3 smooth elliptical patches whose
% total area is motorFraction of the fascicle area (assuming that fibre number is proportional to area) are placed at random,
% without overlap and inside the outline (margin_m from the border). The patches do NOT change the conductivity of the tissue
% (the endoneurial conductivity is an effective value that already includes the axons): they only say where the motor fibres run.
%% Inputs / outputs
%   outline        Nx2 outline of the fascicle [m]
%   motorFraction  fraction (0..1) of the fascicle area covered by motor patches
%   rs             RandStream (reproducibility)
%   margin_m       minimum distance of a patch from the fascicle border [m] (default 10 um)
%   patches        cell array of 24x2 elliptical outlines [m] (fewer than planned if a patch could not be placed in 400 tries)
if nargin<4, margin_m = 10e-6; end
A = polyarea(outline(:,1),outline(:,2)); target = motorFraction*A; patches = {}; % Fascicle area and area to be covered by patches.
if target<=0, return; end
nP = min(3,max(1,round(1+2*rand(rs)))); parts = diff([0;sort(rand(rs,nP-1,1));1])*target; % Number of patches and their areas (random split).
inner = polybuffer(polyshape(outline),-margin_m); th = (0:23)'*2*pi/24; % Allowed region and patch outline resolution.
lo = min(outline,[],1); hi = max(outline,[],1);
for k = 1:nP
    for attempt = 1:400
        c = lo+rand(rs,1,2).*(hi-lo); if ~isinterior(inner,c(1),c(2)), continue; end % Centre inside the fascicle.
        ar = 1+rand(rs); phi = pi*rand(rs); a = sqrt(parts(k)*ar/pi); b = sqrt(parts(k)/(pi*ar)); % Axis ratio 1 to 2, orientation, semi-axes.
        P = [a*cos(th),b*sin(th)]*[cos(phi) sin(phi);-sin(phi) cos(phi)]+c; % Elliptical patch.
        if ~all(isinterior(inner,P(:,1),P(:,2))), continue; end % Stay inside the fascicle.
        clash = false; for j = 1:numel(patches), if overlaps(polyshape(P),polyshape(patches{j})), clash = true; end, end
        if clash, continue; end
        patches{end+1} = P; break %#ok<AGROW> % Accept this patch.
    end
end
end
